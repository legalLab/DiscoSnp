/*
 * disco_haplotypes: C++ implementation of scripts/disco_haplotypes.py
 *
 * Locus-level post-processing of DiscoSnp++ / DiscoSnpRad bubbles: multi-allelic
 * sites, close SNPs and read-backed haplotypes.  Three sub-commands, with the
 * options and the outputs of the python script (the outputs are identical):
 *
 *   augment   before kissreads2: single-SNP "synthetic" bubbles in the sequence
 *             contexts needed by kissreads2
 *   call      after kissreads2 (-phasing, -phasing_sites): sites, genotypes and
 *             haplotypes: <out>.vcf, <out>.tsv, <out>_loci.tsv, <out>_loci.fa,
 *             <out>_alleles.fa
 *   strip     removes the synthetic bubbles from a fasta file
 *
 * See scripts/disco_haplotypes.py for the description of the methods: every
 * function here has the name of the python function it implements.  Every tie
 * is broken as in the python script, so that both give the same files.
 */

#include <algorithm>
#include <array>
#include <climits>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
typedef int64_t i64;
typedef uint64_t u64;

static const uint8_t PAD = 255;
static const int MAX_PATH_LENGTH = 1000;
static const int MAX_FLANK = 1500;
static const int OFFSET_BITS = 13;
static const i64 OFFSET_SHIFT = 1 << (OFFSET_BITS - 1);
static const double ERROR_RATE = 0.01;
static const char NONE = 0;                 // "None" nucleotide
// tolerated mismatches in the overlap of two bubbles: max_mismatches + max_divergence x overlap
static const int STRICT_MISMATCHES = 1;  static const double STRICT_DIVERGENCE = 0.01;
static const int RELAXED_MISMATCHES = 4; static const double RELAXED_DIVERGENCE = 0.02;

static uint8_t CODE[256];
static uint8_t COMP[256];
static const char DECODE[6] = "ACGTN";

static void init_tables() {
    for (int i = 0; i < 256; i++) { CODE[i] = 4; COMP[i] = PAD; }
    const char *s = "ACGTacgt";
    for (int i = 0; i < 8; i++) CODE[(uint8_t) s[i]] = i % 4;
    COMP[0] = 3; COMP[1] = 2; COMP[2] = 1; COMP[3] = 0; COMP[4] = 4;
}

static void log_msg(const string &m) { cerr << m << "\n"; cerr.flush(); }
[[noreturn]] static void fatal(const string &m) { cerr << m << "\n"; exit(1); }

static char complement(char c) {            // str.translate("ACGTNacgtn" -> "TGCANtgcan")
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; case 'N': return 'N';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a'; case 'n': return 'n';
        default: return c;
    }
}
static string revcomp(const string &s) {
    string r(s.rbegin(), s.rend());
    for (auto &c : r) c = complement(c);
    return r;
}
static bool file_exists(const string &p) { struct stat st; return stat(p.c_str(), &st) == 0; }
static string basename_of(const string &p) { auto k = p.find_last_of('/'); return k == string::npos ? p : p.substr(k + 1); }
static bool is_lower(char c) { return c >= 'a' && c <= 'z'; }
static bool py_isupper(char c) { return c >= 'A' && c <= 'Z'; }       // single ASCII characters

// python repr() of a float
static string py_float(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char buf[64];
    auto r = to_chars(buf, buf + sizeof(buf), x, chars_format::scientific);
    string s(buf, r.ptr);
    string sign;
    if (s[0] == '-') { sign = "-"; s = s.substr(1); }
    auto e = s.find('e');
    string mant = s.substr(0, e);
    int exp10 = atoi(s.c_str() + e + 1);
    string digits;
    for (char c : mant) if (c != '.') digits += c;
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    if (x == 0) return sign + "0.0";
    if (exp10 >= -4 && exp10 < 16) {
        string out;
        if (exp10 >= 0) {
            string ip = digits.substr(0, min((size_t) exp10 + 1, digits.size()));
            while ((int) ip.size() < exp10 + 1) ip += '0';
            string fp = (int) digits.size() > exp10 + 1 ? digits.substr(exp10 + 1) : "";
            out = ip + "." + (fp.empty() ? "0" : fp);
        } else {
            out = "0." + string(-exp10 - 1, '0') + digits;
        }
        return sign + out;
    }
    string m = digits.substr(0, 1);
    if (digits.size() > 1) m += "." + digits.substr(1);
    char eb[16];
    snprintf(eb, sizeof(eb), "e%c%02d", exp10 < 0 ? '-' : '+', abs(exp10));
    return sign + m + eb;
}

static const map<string, char> &iupac_table() {
    static map<string, char> t = {{"A", 'A'}, {"C", 'C'}, {"G", 'G'}, {"T", 'T'}, {"AG", 'R'}, {"CT", 'Y'}, {"CG", 'S'},
                                  {"AT", 'W'}, {"GT", 'K'}, {"AC", 'M'}, {"CGT", 'B'}, {"AGT", 'D'}, {"ACT", 'H'},
                                  {"ACG", 'V'}, {"ACGT", 'N'}};
    return t;
}
static char iupac(const vector<char> &alleles) {
    string k;
    for (char c : alleles) if (k.find(c) == string::npos) k += c;
    sort(k.begin(), k.end());
    auto it = iupac_table().find(k);
    return it == iupac_table().end() ? 'N' : it->second;
}

// ------------------------------------------------------------------ input lines
struct LineReader {
    FILE *f; char *buf = nullptr; size_t cap = 0; bool ends_with_newline = true;
    explicit LineReader(const string &path) {
        f = fopen(path.c_str(), "rb");
        if (!f) fatal("ERROR: cannot open " + path);
        setvbuf(f, nullptr, _IOFBF, 1 << 22);
    }
    ~LineReader() { if (f) fclose(f); free(buf); }
    // line without its '\n' (python "line" minus the final newline)
    bool next(string &line) {
        ssize_t n = getline(&buf, &cap, f);
        if (n < 0) return false;
        ends_with_newline = n > 0 && buf[n - 1] == '\n';
        line.assign(buf, ends_with_newline ? n - 1 : n);
        return true;
    }
};

struct Out {
    FILE *f;
    explicit Out(const string &path) {
        f = fopen(path.c_str(), "wb");
        if (!f) fatal("ERROR: cannot write " + path);
        setvbuf(f, nullptr, _IOFBF, 1 << 22);
    }
    ~Out() { fclose(f); }
    Out &operator<<(const string &s) { fwrite(s.data(), 1, s.size(), f); return *this; }
    Out &operator<<(const char *s) { fputs(s, f); return *this; }
    Out &operator<<(char c) { fputc(c, f); return *this; }
    Out &operator<<(i64 v) { char b[32]; auto r = to_chars(b, b + 32, v); fwrite(b, 1, r.ptr - b, f); return *this; }
    Out &operator<<(int v) { return *this << (i64) v; }
    Out &operator<<(size_t v) { return *this << (i64) v; }
};
static string S(i64 v) { char b[32]; auto r = to_chars(b, b + 32, v); return string(b, r.ptr); }

// ------------------------------------------------------------------ headers
// >(SNP|INDEL)_(higher|lower)_path_(\d+)  at the start of the header
static bool header_match(const string &h, bool &snp, bool &higher, i64 &id) {
    size_t p = 1;
    if (h.size() < 2 || h[0] != '>') return false;
    if (h.compare(1, 4, "SNP_") == 0) { snp = true; p = 5; }
    else if (h.compare(1, 6, "INDEL_") == 0) { snp = false; p = 7; }
    else return false;
    if (h.compare(p, 7, "higher_") == 0) { higher = true; p += 7; }
    else if (h.compare(p, 6, "lower_") == 0) { higher = false; p += 6; }
    else return false;
    if (h.compare(p, 5, "path_") != 0) return false;
    p += 5;
    if (p >= h.size() || !isdigit((unsigned char) h[p])) return false;
    id = 0;
    while (p < h.size() && isdigit((unsigned char) h[p])) id = id * 10 + (h[p++] - '0');
    return true;
}
// rank_([^|\s]+), NaN when absent or not a number
static double parse_rank(const string &h) {
    auto k = h.find("rank_");
    if (k == string::npos) return NAN;
    size_t e = k + 5;
    while (e < h.size() && h[e] != '|' && !isspace((unsigned char) h[e])) e++;
    if (e == k + 5) return NAN;
    string tok = h.substr(k + 5, e - k - 5);
    // python float(): digits, sign, dot, exponent, inf, nan (no hexadecimal)
    for (char c : tok) if (!(isdigit((unsigned char) c) || strchr("+-.eEinfatyINFATY", c))) return NAN;
    char *end;
    double v = strtod(tok.c_str(), &end);
    if (*end != 0) return NAN;
    return v;
}
static bool search_pair(const string &h, const char *a, const char *b, int &x, int &y) {
    // a(\d+)\|b(\d+)
    size_t from = 0, la = strlen(a), lb = strlen(b);
    while (true) {
        auto k = h.find(a, from);
        if (k == string::npos) return false;
        size_t p = k + la;
        i64 v1 = 0; size_t d = p;
        while (d < h.size() && isdigit((unsigned char) h[d])) v1 = v1 * 10 + (h[d++] - '0');
        if (d > p && d < h.size() && h[d] == '|' && h.compare(d + 1, lb, b) == 0) {
            size_t q = d + 1 + lb, d2 = q; i64 v2 = 0;
            while (d2 < h.size() && isdigit((unsigned char) h[d2])) v2 = v2 * 10 + (h[d2++] - '0');
            if (d2 > q) { x = (int) v1; y = (int) v2; return true; }
        }
        from = k + 1;
    }
}
static array<int, 4> header_meta(const string &h) {
    array<int, 4> m = {-1, -1, -1, -1};
    int x, y;
    if (search_pair(h, "left_unitig_length_", "right_unitig_length_", x, y)) { m[0] = x; m[1] = y; }
    if (search_pair(h, "left_contig_length_", "right_contig_length_", x, y)) { m[2] = x; m[3] = y; }
    return m;
}
// \|C\d+_(\d+)
static vector<int> header_counts(const string &h) {
    vector<int> v;
    for (size_t k = 0; k + 1 < h.size(); k++) {
        if (h[k] != '|' || h[k + 1] != 'C') continue;
        size_t p = k + 2, d = p;
        while (d < h.size() && isdigit((unsigned char) h[d])) d++;
        if (d == p || d >= h.size() || h[d] != '_') continue;
        size_t q = d + 1, e = q; i64 val = 0;
        while (e < h.size() && isdigit((unsigned char) h[e])) val = val * 10 + (h[e++] - '0');
        if (e == q) continue;
        v.push_back((int) min<i64>(val, 1 << 30));
    }
    return v;
}

// ------------------------------------------------------------------ bubble store
struct Store {
    vector<i64> ids;
    vector<string> paths;                   // [2n] upper-case parts, decoded (ACGT, N)
    vector<u64> foff; vector<uint8_t> fcodes; // [2n] paths with their (capped) extensions, codes
    vector<int> flens, fstart, left, right;
    vector<uint16_t> counts;                // [2n * S]
    vector<double> ranks;
    vector<array<int, 4>> meta;
    int S = 0, n = 0, max_len = 0, max_flen = 0;
    vector<int> by_id;

    const uint8_t *frow(int r) const { return fcodes.data() + foff[r]; }
    int lens(int r) const { return (int) paths[r].size(); }
    void finish() {
        n = (int) ids.size();
        max_len = 0; max_flen = 0;
        for (auto &p : paths) max_len = max(max_len, (int) p.size());
        for (int f : flens) max_flen = max(max_flen, f);
        i64 mx = 0;
        for (i64 i : ids) mx = max(mx, i);
        by_id.assign(n ? mx + 1 : 1, -1);
        for (int i = 0; i < n; i++) by_id[ids[i]] = i;
    }
    int index_of(i64 id) const { return (id < 0 || id >= (i64) by_id.size()) ? -1 : by_id[id]; }
    vector<int> snp_positions(int index) const {
        vector<int> out;
        const string &h = paths[2 * index], &l = paths[2 * index + 1];
        for (size_t p = 0; p < min(h.size(), l.size()); p++) if (h[p] != l[p]) out.push_back((int) p);
        return out;
    }
};

static Store parse_store(const string &file, bool with_counts, const unordered_set<i64> *only_ids, int max_flank) {
    Store st;
    LineReader in(file);
    string line, header, pending_header, pending_path, pending_full;
    bool has_header = false, has_pending = false;
    i64 pending_id = -1;
    int pending_start = 0, pending_left = 0, pending_right = 0;
    int nsamples = -1, n_indel = 0;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '>') { header = line; has_header = true; continue; }
        if (!has_header) continue;
        bool snp, higher; i64 id;
        if (!header_match(header, snp, higher, id)) fatal("ERROR: unexpected header in " + file + ": " + header);
        string this_header = header;
        has_header = false;
        if (!snp) { n_indel += higher; continue; }
        if (only_ids && !only_ids->count(id)) continue;
        // python: the line (with its end of line) without its lower-case letters, \r and \n
        string path;
        for (char c : line) if (!is_lower(c) && c != '\r') path += c;
        for (auto &c : path) c = DECODE[CODE[(uint8_t) c]];
        string seq = line;
        while (!seq.empty() && (seq.back() == '\r' || seq.back() == '\n')) seq.pop_back();
        size_t left = 0;
        while (left < seq.size() && is_lower(seq[left])) left++;
        size_t right = 0;
        while (right < seq.size() && is_lower(seq[seq.size() - 1 - right])) right++;
        if (left == seq.size()) { left = 0; right = 0; }
        int kl = min((int) left, max_flank), kr = min((int) right, max_flank);
        string full = seq.substr(left - kl, seq.size() - right + kr - (left - kl));
        if (higher) {
            pending_id = id; pending_path = path; pending_header = this_header; pending_full = full;
            pending_start = kl; pending_left = (int) left; pending_right = (int) right; has_pending = true;
            continue;
        }
        if (!has_pending || pending_id != id) fatal("ERROR: lower path " + S(id) + " is not preceded by its higher path in " + file);
        if ((int) max(path.size(), pending_path.size()) > MAX_PATH_LENGTH)
            fatal("ERROR: bubble " + S(id) + " has an upper-case path longer than " + S(MAX_PATH_LENGTH) + " nt");
        st.ids.push_back(id);
        auto add = [&](const string &p, const string &f, int start, int l, int r) {
            st.paths.push_back(p);
            st.foff.push_back(st.fcodes.size());
            for (char c : f) st.fcodes.push_back(CODE[(uint8_t) c]);
            st.flens.push_back((int) f.size());
            st.fstart.push_back(start); st.left.push_back(l); st.right.push_back(r);
        };
        add(pending_path, pending_full, pending_start, pending_left, pending_right);
        add(path, full, kl, (int) left, (int) right);
        st.meta.push_back(header_meta(pending_header));
        st.ranks.push_back(parse_rank(pending_header));
        if (with_counts) {
            for (const string *t : {&pending_header, &this_header}) {
                vector<int> v = header_counts(*t);
                for (auto &x : v) x = min(x, 65535);
                if (nsamples < 0) nsamples = (int) v.size();
                else if ((int) v.size() != nsamples)
                    fatal("ERROR: bubble " + S(id) + ": " + S(v.size()) + " read counts, " + S(nsamples) + " expected");
                for (int x : v) st.counts.push_back((uint16_t) x);
            }
        }
        has_pending = false;
    }
    if (has_pending) fatal("ERROR: higher path " + S(pending_id) + " has no lower path in " + file);
    st.finish();
    if (st.max_flen >= (1 << (OFFSET_BITS - 1)))
        fatal("ERROR: paths with their extensions longer than " + S((1 << (OFFSET_BITS - 1)) - 1) + " nt: lower --max_flank");
    if (with_counts) {
        if (st.n && nsamples <= 0) fatal("ERROR: no read counts (C1_, C2_...) in the headers of " + file + ": is this a kissreads2 output?");
        st.S = max(nsamples, 0);
    }
    log_msg("[" + basename_of(file) + "] " + S(st.n) + " SNP bubbles read (" + S(n_indel) + " INDEL bubbles ignored), longest path "
            + S(st.max_len) + " nt (" + S(st.max_flen) + " nt with its extensions)" + (with_counts ? ", " + S(nsamples < 0 ? 0 : nsamples) + " read sets" : ""));
    return st;
}

// this store followed by the bubbles `keep` of `other`
static void extend_store(Store &st, const Store &other, const vector<char> &keep) {
    for (int i = 0; i < other.n; i++) {
        if (!keep[i]) continue;
        st.ids.push_back(other.ids[i]);
        for (int r = 2 * i; r < 2 * i + 2; r++) {
            st.paths.push_back(other.paths[r]);
            st.foff.push_back(st.fcodes.size());
            st.fcodes.insert(st.fcodes.end(), other.frow(r), other.frow(r) + other.flens[r]);
            st.flens.push_back(other.flens[r]); st.fstart.push_back(other.fstart[r]);
            st.left.push_back(other.left[r]); st.right.push_back(other.right[r]);
            for (int s = 0; s < other.S; s++) st.counts.push_back(other.counts[(size_t) r * other.S + s]);
        }
        st.ranks.push_back(other.ranks[i]);
        st.meta.push_back(other.meta[i]);
    }
    st.finish();
}

struct SyntheticMap {
    vector<i64> order;                               // synthetic ids, file order (python dict)
    unordered_map<i64, pair<i64, i64>> parent;       // synthetic id -> (parent id, start)
};
static SyntheticMap read_synthetic_map(const string &file) {
    SyntheticMap m;
    if (file.empty() || !file_exists(file)) return m;
    LineReader in(file);
    string line;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '#') continue;
        istringstream ss(line);
        i64 a, b, c;
        if (!(ss >> a)) continue;
        ss >> b >> c;
        if (!m.parent.count(a)) m.order.push_back(a);
        m.parent[a] = {b, c};
    }
    return m;
}

// ------------------------------------------------------------------ placement edges from the sequences
struct Edges { vector<i64> b1, b2, shift, rel; size_t size() const { return b1.size(); } };

struct View {                                        // placement_view: some bubbles of the store
    int n; const Store *st; vector<int> index;       // bubble i of the view = bubble index[i] of the store
    const uint8_t *row(int r) const { return st->frow(2 * index[r / 2] + r % 2); }
    int flen(int r) const { return st->flens[2 * index[r / 2] + r % 2]; }
    int fstart(int r) const { return st->fstart[2 * index[r / 2] + r % 2]; }
    int ulen(int r) const { return st->lens(2 * index[r / 2] + r % 2); }
};

static inline u64 pack_key(i64 qb, i64 tb, i64 o, int rel) {
    return ((u64) qb << (OFFSET_BITS + 26)) | ((u64) tb << (OFFSET_BITS + 1)) | ((u64) (o + OFFSET_SHIFT) << 1) | (u64) (rel < 0);
}
static inline void unpack_key(u64 k, i64 &qb, i64 &tb, i64 &o, int &rel) {
    qb = (i64) (k >> (OFFSET_BITS + 26));
    tb = (i64) ((k >> (OFFSET_BITS + 1)) & ((1ULL << 25) - 1));
    o = (i64) ((k >> 1) & ((1ULL << OFFSET_BITS) - 1)) - OFFSET_SHIFT;
    rel = (k & 1) ? -1 : 1;
}

// seed_windows: calls emit(code, position) for the selected windows of one (oriented) path
template <class F>
static void seed_windows(const uint8_t *row, int flen, int start, int ulen, bool higher, int k, int sampling, bool query, F emit) {
    if (flen < k) return;
    int last = start + ulen - k;
    uint32_t code = 0; int bad = 0;                     // number of non-ACGT codes in the window
    for (int j = 0; j < flen; j++) {
        uint8_t c = row[j];
        code = (code << 2) | (c & 3);
        if (c >= 4) bad++;
        if (j >= k && row[j - k] >= 4) bad--;
        int pos = j - k + 1;
        if (pos < 0 || bad) continue;
        bool upper = pos >= start && pos <= last;
        if (upper && !query) upper = ((pos - start) % k == 0) || pos == last;
        bool sampled = false;
        if (higher) { uint32_t h = (uint32_t) (code * 2654435761u) >> 16; sampled = h % (uint32_t) sampling == 0; }
        if (upper || sampled) emit(code, pos);
    }
}

static Edges unique_edges(Edges e) {
    vector<size_t> ord(e.size());
    iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](size_t x, size_t y) {
        return tie(e.b1[x], e.b2[x], e.shift[x], e.rel[x]) < tie(e.b1[y], e.b2[y], e.shift[y], e.rel[y]);
    });
    Edges o;
    for (size_t k = 0; k < ord.size(); k++) {
        size_t i = ord[k];
        if (k && e.b1[i] == o.b1.back() && e.b2[i] == o.b2.back() && e.shift[i] == o.shift.back() && e.rel[i] == o.rel.back()) continue;
        o.b1.push_back(e.b1[i]); o.b2.push_back(e.b2[i]); o.shift.push_back(e.shift[i]); o.rel.push_back(e.rel[i]);
    }
    return o;
}

static Edges sequence_edges(const View &v, int k, int max_mismatches, int min_overlap, double max_divergence, int sampling,
                            int max_candidates = 64) {
    Edges empty;
    int n = v.n, max_flen = v.st->max_flen;
    if (n < 2 || max_flen < k) return empty;
    if (n >= (1 << 25)) fatal("ERROR: more than 33 million bubbles are not supported by the edge search");
    // ---- index, sorted by (seed, row, position) as the stable numpy argsort
    struct Seed { uint32_t code; int32_t row; int16_t off; };
    vector<Seed> index;
    for (int r = 0; r < 2 * n; r++)
        seed_windows(v.row(r), v.flen(r), v.fstart(r), v.ulen(r), r % 2 == 0, k, sampling, false,
                     [&](uint32_t code, int pos) { index.push_back({code, r, (int16_t) pos}); });
    sort(index.begin(), index.end(), [](const Seed &a, const Seed &b) {
        return tie(a.code, a.row, a.off) < tie(b.code, b.row, b.off);
    });
    size_t n_repetitive = 0;
    {
        vector<Seed> kept;
        kept.reserve(index.size());
        for (size_t i = 0; i < index.size();) {
            size_t j = i;
            while (j < index.size() && index[j].code == index[i].code) j++;
            if ((int) (j - i) <= max_candidates) kept.insert(kept.end(), index.begin() + i, index.begin() + j);
            else n_repetitive++;
            i = j;
        }
        index.swap(kept);
    }
    vector<uint32_t> codes(index.size());
    for (size_t i = 0; i < index.size(); i++) codes[i] = index[i].code;
    // ---- query, both orientations
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    vector<vector<u64>> keys(nthreads);
    vector<size_t> dedup_at(nthreads, (size_t) 1 << 22);   // bounded buffer per thread, amortised
    vector<u64> hits(nthreads, 0);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int b = 0; b < n; b++) {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        vector<uint8_t> rcrow;
        for (int r = 2 * b; r < 2 * b + 2; r++) {
            int flen = v.flen(r);
            for (int rel : {1, -1}) {
                const uint8_t *row = v.row(r);
                int start = v.fstart(r);
                if (rel == -1) {
                    rcrow.resize(flen);
                    for (int j = 0; j < flen; j++) rcrow[j] = COMP[row[flen - 1 - j]];
                    row = rcrow.data();
                    start = flen - v.fstart(r) - v.ulen(r);
                }
                seed_windows(row, flen, start, v.ulen(r), r % 2 == 0, k, sampling, true, [&](uint32_t code, int pos) {
                    auto lo = lower_bound(codes.begin(), codes.end(), code) - codes.begin();
                    auto hi = upper_bound(codes.begin() + lo, codes.end(), code) - codes.begin();
                    hits[t] += hi - lo;
                    for (auto i = lo; i < hi; i++) {
                        i64 tb = index[i].row / 2;
                        if (tb == b) continue;
                        keys[t].push_back(pack_key(b, tb, pos - index[i].off, rel));
                    }
                    if (keys[t].size() > dedup_at[t]) {
                        // remove the duplicates; next time when the buffer has doubled again
                        sort(keys[t].begin(), keys[t].end());
                        keys[t].erase(unique(keys[t].begin(), keys[t].end()), keys[t].end());
                        dedup_at[t] = max((size_t) 1 << 22, 2 * keys[t].size());
                    }
                });
            }
        }
    }
    u64 n_hits = accumulate(hits.begin(), hits.end(), (u64) 0);
    vector<u64> all;
    for (auto &kv : keys) { all.insert(all.end(), kv.begin(), kv.end()); vector<u64>().swap(kv); }
    sort(all.begin(), all.end());
    all.erase(unique(all.begin(), all.end()), all.end());
    vector<Seed>().swap(index); vector<uint32_t>().swap(codes);
    // ---- verification: best of the 2 x 2 path pairs
    size_t m = all.size();
    vector<int> best_mm(m, 1000000), best_ov(m, 0);
    #pragma omp parallel for schedule(dynamic, 4096)
    for (size_t c = 0; c < m; c++) {
        i64 qb, tb, o; int rel;
        unpack_key(all[c], qb, tb, o, rel);
        i64 a = max<i64>(0, o), bb = max<i64>(0, -o);
        i64 width = max_flen - max(a, bb);
        if (width <= 0) continue;
        int gmm = 1000000, gov = 0;
        for (int qi = 0; qi < 2; qi++) {
            int rq = 2 * (int) qb + qi, lq = v.flen(rq);
            const uint8_t *Q = v.row(rq);
            for (int ti = 0; ti < 2; ti++) {
                int rt = 2 * (int) tb + ti, lt = v.flen(rt);
                const uint8_t *T = v.row(rt);
                i64 ov = min<i64>(width, min<i64>(lq - a, lt - bb));
                if (ov < 0) ov = 0;
                int mm = 0;
                if (rel == 1) { for (i64 j = 0; j < ov; j++) mm += Q[a + j] != T[bb + j]; }
                else { for (i64 j = 0; j < ov; j++) mm += COMP[Q[lq - 1 - (a + j)]] != T[bb + j]; }
                if (mm < gmm || (mm == gmm && ov > gov)) { gmm = mm; gov = (int) ov; }
            }
        }
        best_mm[c] = gmm; best_ov[c] = gov;
    }
    struct Cand { i64 qb, tb, shift; int rel, mm, ov; };
    vector<Cand> acc;
    for (size_t c = 0; c < m; c++) {
        int mm = best_mm[c], ov = best_ov[c];
        if (!(ov >= min_overlap && (i64) mm * 10 <= ov && mm <= max_mismatches + floor(max_divergence * ov))) continue;
        i64 qb, tb, o; int rel;
        unpack_key(all[c], qb, tb, o, rel);
        i64 fq = v.fstart(2 * (int) qb), ft = v.fstart(2 * (int) tb);
        i64 shift = rel == 1 ? o + ft - fq : (i64) v.flen(2 * (int) qb) - 1 - o - ft - fq;
        acc.push_back({qb, tb, shift, rel, mm, ov});
    }
    if (acc.empty()) return empty;
    // one placement per bubble pair: lexsort((-ov, mm, tb, qb)), stable
    stable_sort(acc.begin(), acc.end(), [](const Cand &x, const Cand &y) {
        return make_tuple(x.qb, x.tb, x.mm, -x.ov) < make_tuple(y.qb, y.tb, y.mm, -y.ov);
    });
    Edges e;
    for (size_t i = 0; i < acc.size(); i++) {
        if (i && acc[i].qb == acc[i - 1].qb && acc[i].tb == acc[i - 1].tb) continue;
        const Cand &x = acc[i];
        bool swap_ = x.qb > x.tb;
        e.b1.push_back(swap_ ? x.tb : x.qb); e.b2.push_back(swap_ ? x.qb : x.tb);
        e.shift.push_back(swap_ ? -x.rel * x.shift : x.shift); e.rel.push_back(x.rel);
    }
    e = unique_edges(e);
    log_msg("[edges] " + S(e.size()) + " sequence overlaps between bubbles (" + S(n_hits) + " seed hits, " + S(n_repetitive)
            + " repetitive seeds ignored)");
    return e;
}

// bubbles overlapping more than max_overlaps other bubbles are repeats: their edges are removed
static Edges remove_repeats(const Edges &e, int n, int max_overlaps, vector<i64> &degree, vector<char> &repeat) {
    degree.assign(n, 0);
    for (size_t i = 0; i < e.size(); i++) { degree[e.b1[i]]++; degree[e.b2[i]]++; }
    repeat.assign(n, 0);
    i64 n_repeat = 0;
    if (max_overlaps > 0) for (int i = 0; i < n; i++) if (degree[i] > max_overlaps) { repeat[i] = 1; n_repeat++; }
    Edges out;
    for (size_t i = 0; i < e.size(); i++) {
        if (repeat[e.b1[i]] || repeat[e.b2[i]]) continue;
        out.b1.push_back(e.b1[i]); out.b2.push_back(e.b2[i]); out.shift.push_back(e.shift[i]); out.rel.push_back(e.rel[i]);
    }
    if (n_repeat)
        log_msg("[repeats] " + S(n_repeat) + " bubbles overlap more than " + S(max_overlaps) + " other bubbles: repeats, kept as isolated bubbles ("
                + S(e.size() - out.size()) + " placement edges)");
    return out;
}

// ------------------------------------------------------------------ facts
struct Token { int sign; i64 id; int path; i64 gap; };

static bool parse_int(const string &s, size_t &p, i64 &v, bool allow_sign) {
    bool neg = false;
    if (allow_sign && p < s.size() && s[p] == '-') { neg = true; p++; }
    size_t d = p; v = 0;
    while (p < s.size() && isdigit((unsigned char) s[p])) v = v * 10 + (s[p++] - '0');
    if (p == d) return false;
    if (neg) v = -v;
    return true;
}
// '-1h_0;-3h_-51;' -> tokens; false when a token does not match (python: None)
static bool parse_fact_part(const string &text, vector<Token> &tokens) {
    tokens.clear();
    size_t b = 0, e = text.size();
    while (b < e && isspace((unsigned char) text[b])) b++;
    while (e > b && isspace((unsigned char) text[e - 1])) e--;
    size_t p = b;
    while (p <= e) {
        size_t q = text.find(';', p);
        if (q == string::npos || q > e) q = e;
        if (q > p) {
            string tok = text.substr(p, q - p);
            size_t x = 0; Token t{1, 0, 0, 0};
            if (x < tok.size() && tok[x] == '-') { t.sign = -1; x++; }
            i64 id;
            if (!parse_int(tok, x, id, false)) return false;
            t.id = id;
            if (x >= tok.size() || (tok[x] != 'h' && tok[x] != 'l')) return false;
            t.path = tok[x] == 'h' ? 0 : 1; x++;
            if (x >= tok.size() || tok[x] != '_') return false;
            x++;
            if (!parse_int(tok, x, t.gap, true)) return false;
            if (x != tok.size()) return false;
            tokens.push_back(t);
        }
        p = q + 1;
    }
    return true;
}
// text.split(): whitespace separated parts
static void split_ws(const string &s, vector<string> &out) {
    out.clear();
    size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && isspace((unsigned char) s[i])) i++;
        size_t j = i;
        while (j < s.size() && !isspace((unsigned char) s[j])) j++;
        if (j > i) out.push_back(s.substr(i, j - i));
        i = j;
    }
}
// parts of a classic fact line (the parts that parse into at least one token)
static void fact_parts(const string &text, vector<vector<Token>> &parts) {
    vector<string> words;
    split_ws(text, words);
    parts.clear();
    vector<Token> t;
    for (auto &w : words) if (parse_fact_part(w, t) && !t.empty()) parts.push_back(t);
}
static bool fact_line(const string &line, string &text, i64 &support) {
    if (!line.empty() && line[0] == '#') return false;
    auto k = line.rfind("=>");
    if (k == string::npos) return false;
    text = line.substr(0, k);
    string s = line.substr(k + 2);
    char *end;
    support = strtoll(s.c_str(), &end, 10);
    return true;
}
static int sample_of_fact_file(const string &f) {
    string b = basename_of(f);
    auto k = b.find("read_set_id_");
    if (k == string::npos || k + 12 >= b.size() || !isdigit((unsigned char) b[k + 12]))
        fatal("ERROR: cannot find the read set index in the name of " + f);
    return atoi(b.c_str() + k + 12) - 1;
}
static string read_set_name(const string &f) {
    LineReader in(f);
    string line;
    if (!in.next(line)) return ".";
    if (line.empty() || line[0] != '#') return ".";
    string s = line.substr(1);
    size_t b = 0, e = s.size();
    while (b < e && isspace((unsigned char) s[b])) b++;
    while (e > b && isspace((unsigned char) s[e - 1])) e--;
    return s.substr(b, e - b);
}

static Edges fact_edges(const vector<string> &fact_files, const Store &st) {
    set<tuple<i64, i64, i64, i64>> keys;
    unordered_set<string> seen;
    i64 n_unknown = 0;
    string line, text;
    vector<vector<Token>> parts;
    for (auto &file : fact_files) {
        LineReader in(file);
        i64 support;
        while (in.next(line)) {
            line += '\n';
            if (!fact_line(line, text, support)) continue;
            if (!seen.insert(text).second) continue;
            fact_parts(text, parts);
            for (auto &part : parts) {
                bool has_prev = false;
                i64 p_start = 0, p_len = 0, p_zero = 0; int p_sign = 0, p_index = 0;
                i64 start = 0;
                for (auto &t : part) {
                    int index = st.index_of(t.id);
                    if (index < 0) { n_unknown++; has_prev = false; continue; }
                    i64 length = st.lens(2 * index + t.path);
                    if (has_prev) start = p_start + p_len + t.gap;
                    i64 zero = t.sign == 1 ? start : start + length - 1;
                    if (has_prev) {
                        i64 shift = p_sign * (zero - p_zero), rel = p_sign * t.sign;
                        i64 b1 = p_index < index ? p_index : index, b2 = p_index < index ? index : p_index;
                        if (p_index > index) shift = -rel * shift;
                        keys.insert({b1, b2, shift, rel});
                    }
                    has_prev = true; p_start = start; p_len = length; p_sign = t.sign; p_index = index; p_zero = zero;
                }
            }
        }
    }
    if (n_unknown) log_msg("[facts] " + S(n_unknown) + " fact tokens name bubbles absent from the fasta files (ignored)");
    Edges e;
    for (auto &k : keys) { e.b1.push_back(get<0>(k)); e.b2.push_back(get<1>(k)); e.shift.push_back(get<2>(k)); e.rel.push_back(get<3>(k)); }
    return e;
}

// ------------------------------------------------------------------ loci
struct Snp { int pos; char high, low; };

struct Bubble {
    i64 id; int index;
    string paths[2];
    vector<Snp> snps;
    string headers[2], sequences[2];
    bool has_headers = false;
    i64 parent = -1; bool synthetic = false;
    i64 anchor = 0; int orient = 1; int locus = -1;
    int left = 0, right = 0;
    string full;
    array<int, 4> meta;
    double rank;
    i64 coordinate(i64 position) const { return anchor + orient * position; }
    i64 position(i64 coordinate) const { return orient * (coordinate - anchor); }
    char nucleotide(int path_index, i64 coord) const {
        i64 p = position(coord);
        const string &path = paths[path_index];
        if (p < 0 || p >= (i64) path.size()) return NONE;
        return orient == 1 ? path[p] : complement(path[p]);
    }
    pair<i64, i64> extent() const { return {-left, (i64) max(paths[0].size(), paths[1].size()) + right - 1}; }
    bool has_rank() const { return !std::isnan(rank); }
};

struct Carrier { int bubble; int path; };
struct SiteCarriers {                                 // nucleotide -> carriers, in order of first appearance
    vector<pair<char, vector<Carrier>>> v;
    vector<Carrier> &operator[](char n) {
        for (auto &x : v) if (x.first == n) return x.second;
        v.push_back({n, {}});
        return v.back().second;
    }
    const vector<Carrier> *get(char n) const {
        for (auto &x : v) if (x.first == n) return &x.second;
        return nullptr;
    }
};
struct Locus {
    string name;
    vector<int> bubbles;                              // slots, sorted by store index
    i64 conflicts = 0, length = 0;
    map<i64, SiteCarriers> sites;
};

struct Materialised {
    vector<Bubble> bubbles;                           // slots
    vector<int> slot_of;                              // store index -> slot, -1
    unordered_map<i64, int> slot_of_id;
};

static Bubble make_bubble(const Store &st, int index) {
    Bubble b;
    b.id = st.ids[index]; b.index = index;
    b.paths[0] = st.paths[2 * index]; b.paths[1] = st.paths[2 * index + 1];
    const string &h = b.paths[0], &l = b.paths[1];
    for (size_t p = 0; p < min(h.size(), l.size()); p++) if (h[p] != l[p]) b.snps.push_back({(int) p, h[p], l[p]});
    b.left = st.left[2 * index]; b.right = st.right[2 * index];
    b.meta = st.meta[index]; b.rank = st.ranks[index];
    return b;
}

static vector<Locus> build_loci(Materialised &M, const Edges &edges) {
    int nb = (int) M.bubbles.size();
    vector<vector<tuple<int, i64, i64>>> graph(nb);
    for (size_t i = 0; i < edges.size(); i++) {
        int a = M.slot_of[edges.b1[i]], b = M.slot_of[edges.b2[i]];
        if (a < 0 || b < 0) continue;
        graph[a].push_back({b, edges.shift[i], edges.rel[i]});
        graph[b].push_back({a, -edges.rel[i] * edges.shift[i], edges.rel[i]});
    }
    vector<int> order(nb);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y) { return M.bubbles[x].index < M.bubbles[y].index; });
    vector<Locus> loci;
    for (int seed : order) {
        Bubble &sb = M.bubbles[seed];
        if (sb.locus >= 0) continue;
        int li = (int) loci.size();
        loci.emplace_back();
        Locus &locus = loci.back();
        sb.locus = li; sb.anchor = 0; sb.orient = 1;
        vector<int> stack = {seed};
        while (!stack.empty()) {
            int cur = stack.back(); stack.pop_back();
            locus.bubbles.push_back(cur);
            const Bubble &c = M.bubbles[cur];
            for (auto &[nbr, shift, relative] : graph[cur]) {
                Bubble &nb_ = M.bubbles[nbr];
                i64 anchor = c.anchor + c.orient * shift;
                int orient = (int) (c.orient * relative);
                if (nb_.locus < 0) { nb_.locus = li; nb_.anchor = anchor; nb_.orient = orient; stack.push_back(nbr); }
                else if (nb_.anchor != anchor || nb_.orient != orient) locus.conflicts++;
            }
        }
    }
    for (auto &locus : loci) {
        locus.conflicts /= 2;
        sort(locus.bubbles.begin(), locus.bubbles.end(), [&](int x, int y) { return M.bubbles[x].index < M.bubbles[y].index; });
        i64 lowest = LLONG_MAX;
        for (int s : locus.bubbles) {
            auto ex = M.bubbles[s].extent();
            lowest = min({lowest, M.bubbles[s].coordinate(ex.first), M.bubbles[s].coordinate(ex.second)});
        }
        i64 highest = 0;
        for (int s : locus.bubbles) {
            Bubble &b = M.bubbles[s];
            b.anchor -= lowest;
            auto ex = b.extent();
            highest = max({highest, b.coordinate(ex.first), b.coordinate(ex.second)});
        }
        locus.length = highest + 1;
        for (int s : locus.bubbles) {
            Bubble &b = M.bubbles[s];
            for (auto &snp : b.snps) {
                for (int pi = 0; pi < 2; pi++) {
                    char nt = pi == 0 ? snp.high : snp.low;
                    if (b.orient == -1) nt = complement(nt);
                    locus.sites[b.coordinate(snp.pos)][nt].push_back({s, pi});
                }
            }
        }
    }
    return loci;
}

// read_link_agrees: the extended higher paths of a and b, placed as the read link says, agree as two bubbles
// placed by sequence must (same acceptance rule as sequence_edges)
static bool read_link_agrees(const Store &st, i64 a, i64 b, i64 shift, i64 rel, int min_overlap, int max_mismatches,
                             double max_divergence) {
    i64 la = st.flens[2 * a], sa = st.fstart[2 * a], lb = st.flens[2 * b], sb = st.fstart[2 * b];
    const uint8_t *fa = st.frow(2 * a), *fb = st.frow(2 * b);
    i64 base = shift - rel * sb + sa, j0, j1;
    if (rel == 1) { j0 = max<i64>(0, -base); j1 = min<i64>(lb, la - base); }
    else { j0 = max<i64>(0, base - la + 1); j1 = min<i64>(lb, base + 1); }
    if (j1 <= j0) return false;
    i64 overlap = 0, mismatches = 0;
    for (i64 j = j0; j < j1; j++) {
        uint8_t xb = rel == 1 ? fb[j] : COMP[fb[j]], xa = fa[base + rel * j];
        if (xa < 4 && xb < 4) { overlap++; mismatches += xa != xb; }
    }
    return overlap >= min_overlap && mismatches * 10 <= overlap && mismatches <= max_mismatches + floor(max_divergence * overlap);
}

// the placement edges that build the loci (see the python limit_locus_edges): read links whose bubbles disagree
// in sequence are removed; with max_locus_length > 0, a union-find keeping the placement of every bubble in the
// frame of its root and the extent of every component takes the synthetic -> parent edges, then the sequence
// overlaps, then the read links, rejecting an edge that makes a longer locus or contradicts a placement
enum { PARENT_EDGE = 0, SEQUENCE_EDGE = 1, READ_EDGE = 2 };
struct EdgeRemovals { i64 read_disagree = 0, sequence_long = 0, read_long = 0, sequence_contradiction = 0, read_contradiction = 0; };

static vector<char> limit_locus_edges(const Edges &edges, const vector<char> &kind, const Store &st, int max_locus_length,
                                      int min_overlap, int max_mismatches, double max_divergence, EdgeRemovals &removed,
                                      vector<i64> &contradictions) {
    vector<char> keep(edges.size(), 1);
    removed = EdgeRemovals();
    contradictions.clear();
    for (size_t i = 0; i < edges.size(); i++)
        if (kind[i] == READ_EDGE && !read_link_agrees(st, edges.b1[i], edges.b2[i], edges.shift[i], edges.rel[i], min_overlap,
                                                      max_mismatches, max_divergence)) {
            keep[i] = 0; removed.read_disagree++;
        }
    if (max_locus_length <= 0) return keep;
    vector<i64> parent(st.n, -1), offset(st.n, 0), low(st.n, 0), high(st.n, 0);
    vector<int> orient(st.n, 1);
    vector<i64> path;
    auto find_root = [&](i64 x) -> i64 {
        if (parent[x] < 0) {
            parent[x] = x; offset[x] = 0; orient[x] = 1;
            low[x] = -(i64) st.left[2 * x];
            high[x] = (i64) max(st.lens(2 * x), st.lens(2 * x + 1)) + st.right[2 * x] - 1;
            return x;
        }
        path.clear();
        while (parent[x] != x) { path.push_back(x); x = parent[x]; }
        for (auto it = path.rbegin(); it != path.rend(); ++it) {     // nearest to the root first
            i64 node = *it, up = parent[node];
            if (up != x) {
                offset[node] = offset[up] + orient[up] * offset[node];
                orient[node] = orient[up] * orient[node];
                parent[node] = x;
            }
        }
        return x;
    };
    // 0: accepted, 1: longer locus, 2: contradicts the placement
    auto unite = [&](i64 a, i64 b, i64 shift, i64 rel, bool check) -> int {
        i64 ra = find_root(a), rb = find_root(b);
        i64 oa = a == ra ? 0 : offset[a], ob = b == rb ? 0 : offset[b];
        int sa = a == ra ? 1 : orient[a], sb = b == rb ? 1 : orient[b];
        if (ra == rb) return check && (oa + sa * shift != ob || sa * rel != sb) ? 2 : 0;
        int o = (int) (sa * rel * sb);
        i64 t = oa + sa * shift - o * ob;
        i64 first = t + o * low[rb], last = t + o * high[rb];
        if (first > last) swap(first, last);
        i64 merged_low = min(low[ra], first), merged_high = max(high[ra], last);
        if (check && merged_high - merged_low + 1 > max_locus_length) return 1;
        parent[rb] = ra; offset[rb] = t; orient[rb] = o;
        low[ra] = merged_low; high[ra] = merged_high;
        return 0;
    };
    for (int current : {PARENT_EDGE, SEQUENCE_EDGE, READ_EDGE})
        for (size_t i = 0; i < edges.size(); i++) {
            if (!keep[i] || kind[i] != current) continue;
            int why = unite(edges.b1[i], edges.b2[i], edges.shift[i], edges.rel[i], current != PARENT_EDGE);
            if (!why) continue;
            keep[i] = 0;
            if (why == 1) (current == SEQUENCE_EDGE ? removed.sequence_long : removed.read_long)++;
            else {
                (current == SEQUENCE_EDGE ? removed.sequence_contradiction : removed.read_contradiction)++;
                contradictions.push_back(edges.b1[i]);
            }
        }
    return keep;
}

// materialise the bubbles of the components (2..max_locus_bubbles bubbles), build their loci
// synthetic: flags of the synthetic bubbles (store indices), not counted in the size of a component
static vector<Locus> multi_bubble_loci(const Store &st, const Edges &edges, int max_locus_bubbles, bool lone_multi_snp,
                                       Materialised &M, const vector<char> *synthetic = nullptr,
                                       const vector<char> *excluded = nullptr, vector<char> *oversized_out = nullptr) {
    vector<int> parent(st.n, -1);
    function<int(int)> find = [&](int x) {
        int r = x;
        while (parent[r] != r) r = parent[r];
        while (parent[x] != r) { int nx = parent[x]; parent[x] = r; x = nx; }
        return r;
    };
    for (size_t i = 0; i < edges.size(); i++) {
        int a = (int) edges.b1[i], b = (int) edges.b2[i];
        if (parent[a] < 0) parent[a] = a;
        if (parent[b] < 0) parent[b] = b;
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
    }
    vector<int> size(st.n, 0);
    for (int i = 0; i < st.n; i++) if (parent[i] >= 0 && !(synthetic && (*synthetic)[i])) size[find(i)]++;
    M.slot_of.assign(st.n, -1);
    auto add = [&](int index) {
        if (M.slot_of[index] >= 0) return;
        M.slot_of[index] = (int) M.bubbles.size();
        M.bubbles.push_back(make_bubble(st, index));
        M.slot_of_id[st.ids[index]] = M.slot_of[index];
    };
    i64 n_oversized = 0, n_components = 0;
    if (lone_multi_snp)
        for (int i = 0; i < st.n; i++)
            if (parent[i] < 0 && st.snp_positions(i).size() > 1 && !(excluded && (*excluded)[i])) add(i);
    for (int i = 0; i < st.n; i++) {
        if (parent[i] < 0) continue;
        int r = find(i);
        if (size[r] > max_locus_bubbles) {
            n_oversized++;
            if (r == i) n_components++;
            if (oversized_out) (*oversized_out)[i] = 1;
            continue;
        }
        add(i);
    }
    vector<Locus> loci = build_loci(M, edges);
    if (n_oversized)
        log_msg("[loci] " + S(n_oversized) + " bubbles belong to " + S(n_components) + " components larger than --max_locus_bubbles ("
                + S(max_locus_bubbles) + "): probably repeats, treated as isolated bubbles");
    return loci;
}

// ------------------------------------------------------------------ options
struct Args {
    string command, input, output, map, coherent, uncoherent, out;
    vector<string> phased, sites;
    int max_contexts = 32, seed_window = 31, min_depth = 3, min_support = 2, min_sites = 1, max_locus_length = 0;
    int seed_size = 16, max_mismatches = -1, min_overlap = 25, max_flank = 1000, seed_sampling = 8, max_locus_bubbles = 500;
    double max_divergence = -1;                       // -1: default (strict, or relaxed with --relaxed_placement)
    bool relaxed_placement = false;
    int max_overlaps = 50;                           // see the python help
};

// ------------------------------------------------------------------ augment
static pair<string, string> canonical_pair(const string &a, const string &b) {
    string x = min(a, revcomp(a)), y = min(b, revcomp(b));
    return x < y ? make_pair(x, y) : make_pair(y, x);
}

struct SubBubble { i64 start; string label; string seq[2]; bool capped; };

static vector<SubBubble> sub_bubble(const Bubble &bubble, int snp_index, const Locus &locus, int max_contexts, int seed_window) {
    vector<SubBubble> out;
    i64 position = bubble.snps[snp_index].pos;
    i64 start = position - bubble.snps[0].pos;
    i64 site = bubble.coordinate(position);
    string lefts[2], rights[2], templates[2];
    for (int pi = 0; pi < 2; pi++) {
        const string &path = bubble.paths[pi];
        i64 stop = position + ((i64) path.size() - 1 - bubble.snps.back().pos);
        const string &full = bubble.sequences[pi];
        i64 first_upper = -1, last_upper = -1;
        for (size_t i = 0; i < full.size(); i++) if (py_isupper(full[i])) { if (first_upper < 0) first_upper = i; last_upper = i; }
        lefts[pi] = start == 0 ? full.substr(0, first_upper) : "";
        rights[pi] = stop == (i64) path.size() - 1 ? full.substr(last_upper + 1) : "";
        templates[pi] = path.substr(start, stop - start + 1);
    }
    i64 shortest = min(templates[0].size(), templates[1].size());
    vector<i64> inside;
    for (auto &kv : locus.sites) {
        i64 c = kv.first;
        i64 rp = bubble.position(c) - start;
        if (c != site && rp >= 0 && rp < shortest) inside.push_back(c);
    }
    stable_sort(inside.begin(), inside.end(), [&](i64 x, i64 y) { return llabs(x - site) < llabs(y - site); });
    vector<i64> blocked;
    for (i64 c : inside) blocked.push_back(bubble.position(c) - start);
    sort(blocked.begin(), blocked.end());
    i64 max_gap = 0, prev = -1;
    for (i64 b : blocked) { max_gap = max(max_gap, b - prev - 1); prev = b; }
    max_gap = max(max_gap, shortest - prev - 1);
    auto build = [&](const vector<i64> &used, const vector<char> &context) {
        SubBubble sb;
        sb.start = start;
        for (int pi = 0; pi < 2; pi++) {
            string path = templates[pi];
            for (size_t u = 0; u < used.size(); u++) {
                char nt = context[u];
                if (nt == NONE) continue;
                if (bubble.orient == -1) nt = complement(nt);
                path[bubble.position(used[u]) - start] = nt;
            }
            path[position - start] = pi == 0 ? bubble.snps[snp_index].high : bubble.snps[snp_index].low;
            sb.seq[pi] = lefts[pi] + path + rights[pi];
        }
        return sb;
    };
    if (max_gap >= seed_window) {
        vector<char> ctx;
        for (i64 c : inside) ctx.push_back(bubble.nucleotide(0, c));
        SubBubble sb = build(inside, ctx);
        sb.label = "."; sb.capped = false;
        out.push_back(sb);
        return out;
    }
    vector<i64> used;
    i64 n_contexts = 1;
    bool capped = false;
    for (i64 c : inside) {
        i64 na = (i64) locus.sites.at(c).v.size();
        if (n_contexts * na > max_contexts) { capped = true; break; }
        n_contexts *= na;
        used.push_back(c);
    }
    vector<vector<char>> choices;
    for (i64 c : used) {
        vector<char> nts;
        for (auto &x : locus.sites.at(c).v) nts.push_back(x.first);
        sort(nts.begin(), nts.end());
        choices.push_back(nts);
    }
    vector<size_t> idx(used.size(), 0);
    while (true) {                                    // itertools.product, last position varying fastest
        vector<char> ctx;
        for (size_t u = 0; u < used.size(); u++) ctx.push_back(choices[u][idx[u]]);
        SubBubble sb = build(used, ctx);
        string label;
        for (size_t u = 0; u < used.size(); u++) {
            if (u) label += ",";
            label += S(used[u] + 1) + ":" + string(1, ctx[u]);
        }
        sb.label = label; sb.capped = capped;
        out.push_back(sb);
        int u = (int) used.size() - 1;
        while (u >= 0 && ++idx[u] == choices[u].size()) { idx[u] = 0; u--; }
        if (u < 0) break;
    }
    return out;
}

static void augment(const Args &args) {
    Store st = parse_store(args.input, false, nullptr, args.max_flank);
    View v{st.n, &st, {}};
    v.index.resize(st.n);
    iota(v.index.begin(), v.index.end(), 0);
    Edges edges = sequence_edges(v, args.seed_size, args.max_mismatches, args.min_overlap, args.max_divergence, args.seed_sampling);
    vector<i64> degree; vector<char> repeat;
    edges = remove_repeats(edges, st.n, args.max_overlaps, degree, repeat);
    Materialised M;
    vector<Locus> all = multi_bubble_loci(st, edges, args.max_locus_bubbles, true, M, nullptr, &repeat);
    vector<Locus> loci;
    for (auto &l : all) if (l.sites.size() > 1) loci.push_back(l);
    for (size_t i = 0; i < loci.size(); i++) loci[i].name = "locus_" + S(i + 1);
    unordered_map<i64, int> wanted;                   // bubble id -> slot
    for (auto &l : loci) for (int s : l.bubbles) wanted[M.bubbles[s].id] = s;
    i64 next_id = 0;
    Out fasta(args.output);
    {
        LineReader in(args.input);
        string line, header;
        bool has_header = false;
        while (in.next(line)) {
            fasta << line;
            if (in.ends_with_newline) fasta << '\n';
            if (!line.empty() && line[0] == '>') { header = line; has_header = true; continue; }
            if (!has_header) continue;
            bool snp, higher; i64 id;
            bool ok = header_match(header, snp, higher, id);
            if (ok) {
                next_id = max(next_id, id + 1);
                auto it = wanted.find(id);
                if (it != wanted.end()) {
                    Bubble &b = M.bubbles[it->second];
                    int pi = higher ? 0 : 1;
                    b.has_headers = true;
                    b.headers[pi] = header;
                    b.sequences[pi] = line;
                }
            }
            has_header = false;
        }
    }
    i64 n_synthetic = 0, n_capped = 0;
    Out mapping(args.map);
    mapping << "#synthetic_id\tparent_id\tstart_in_parent\tlocus\tsite\tcontext\n";
    for (auto &locus : loci) {
        set<pair<string, string>> existing;
        for (int s : locus.bubbles) existing.insert(canonical_pair(M.bubbles[s].paths[0], M.bubbles[s].paths[1]));
        for (int s : locus.bubbles) {
            const Bubble &b = M.bubbles[s];
            for (size_t si = 0; si < b.snps.size(); si++) {
                bool capped = false;
                for (auto &sb : sub_bubble(b, (int) si, locus, args.max_contexts, args.seed_window)) {
                    capped = sb.capped;
                    string up[2];
                    for (int pi = 0; pi < 2; pi++) for (char c : sb.seq[pi]) if (py_isupper(c)) up[pi] += c;
                    auto pair_ = canonical_pair(up[0], up[1]);
                    if (existing.count(pair_)) continue;
                    existing.insert(pair_);
                    for (int pi = 0; pi < 2; pi++) {
                        const string &h = b.headers[pi];
                        // re.sub(r"^>SNP_(higher|lower)_path_\d+", ">SNP_<level>_path_<next_id>")
                        size_t p = h.find("_path_");
                        size_t e = p + 6;
                        while (e < h.size() && isdigit((unsigned char) h[e])) e++;
                        string nh = string(">SNP_") + (pi == 0 ? "higher" : "lower") + "_path_" + S(next_id) + h.substr(e);
                        fasta << nh << "\n" << sb.seq[pi] << "\n";
                    }
                    mapping << next_id << "\t" << b.id << "\t" << sb.start << "\t" << locus.name << "\t"
                            << (b.coordinate(b.snps[si].pos) + 1) << "\t" << (sb.label.empty() ? string(".") : sb.label) << "\n";
                    next_id++; n_synthetic++;
                }
                n_capped += capped;
            }
        }
    }
    log_msg("[augment] " + S(loci.size()) + " loci with several sites (" + S(wanted.size()) + " bubbles)\n[augment] "
            + S(n_synthetic) + " synthetic single-SNP context bubbles written (" + S(n_capped) + " SNPs limited by --max_contexts)");
}

// ------------------------------------------------------------------ call: genotyping
struct Model { vector<pair<int, int>> genotypes; vector<vector<double>> W; };
static Model build_model(int A) {
    Model m;
    for (int b = 0; b < A; b++) for (int a = 0; a <= b; a++) m.genotypes.push_back({a, b});
    for (auto &g : m.genotypes) {
        vector<double> w(A);
        for (int i = 0; i < A; i++) {
            double pa = i == g.first ? 1 - ERROR_RATE : ERROR_RATE / 3;
            double pb = i == g.second ? 1 - ERROR_RATE : ERROR_RATE / 3;
            w[i] = log10((pa + pb) / 2);
        }
        m.W.push_back(w);
    }
    return m;
}
// built once, before any thread starts (read-only afterwards: no data race)
static const Model MODELS[5] = {Model(), Model(), build_model(2), build_model(3), build_model(4)};
static const Model &model(int A) { return MODELS[A]; }
// one site of one sample: gt (-1 when missing), PL, GQ
static void call_genotype(const uint16_t *depths, int A, int min_depth, int gt[2], vector<uint16_t> &PL, int &GQ) {
    const Model &m = model(A);
    size_t G = m.genotypes.size();
    vector<double> LL(G, 0.0);
    for (int i = 0; i < A; i++) {
        double d = depths[i];
        for (size_t g = 0; g < G; g++) LL[g] = LL[g] + d * m.W[g][i];
    }
    double best = *max_element(LL.begin(), LL.end());
    PL.assign(G, 0);
    for (size_t g = 0; g < G; g++) {
        double v = nearbyint(-10 * (LL[g] - best));
        PL[g] = (uint16_t) min(v, 9999.0);
    }
    size_t which = min_element(PL.begin(), PL.end()) - PL.begin();
    vector<uint16_t> sorted_pl = PL;
    sort(sorted_pl.begin(), sorted_pl.end());
    GQ = min<int>(sorted_pl[1], 99);
    gt[0] = m.genotypes[which].first; gt[1] = m.genotypes[which].second;
    i64 total = 0;
    for (int i = 0; i < A; i++) total += depths[i];
    if (total < min_depth) { gt[0] = gt[1] = -1; GQ = 0; }
}
static string join_ints(const vector<i64> &v) {
    string s;
    for (size_t i = 0; i < v.size(); i++) { if (i) s += ","; s += S(v[i]); }
    return s;
}
static string format_field(const int gt[2], bool phased, i64 ps, const vector<i64> &depths, int gq, const vector<i64> &pl) {
    i64 sum = 0;
    for (i64 d : depths) sum += d;
    string ad = join_ints(depths);
    if (gt[0] < 0) return "./.:.:" + S(sum) + ":" + ad + ":.:.";
    if (phased) return S(gt[0]) + "|" + S(gt[1]) + ":" + S(ps) + ":" + S(sum) + ":" + ad + ":" + S(gq) + ":" + join_ints(pl);
    return S(min(gt[0], gt[1])) + "/" + S(max(gt[0], gt[1])) + ":.:" + S(sum) + ":" + ad + ":" + S(gq) + ":" + join_ints(pl);
}

// ------------------------------------------------------------------ call: phasing
struct ParityUnionFind {
    map<i64, i64> parent; map<i64, int> parity;
    explicit ParityUnionFind(const vector<i64> &items) { for (i64 i : items) { parent[i] = i; parity[i] = 0; } }
    pair<i64, int> find(i64 item) {
        vector<i64> path;
        while (parent[item] != item) { path.push_back(item); item = parent[item]; }
        int p = 0;
        for (auto it = path.rbegin(); it != path.rend(); ++it) { p ^= parity[*it]; parent[*it] = item; parity[*it] = p; }
        return {item, p};
    }
    void unite(i64 a, i64 b, int flipped) {
        auto [ra, pa] = find(a);
        auto [rb, pb] = find(b);
        if (ra == rb) return;
        parent[rb] = ra;
        parity[rb] = pa ^ pb ^ flipped;
    }
};

// votes of one locus: (c1, c2, n1, n2) -> support
typedef map<tuple<i64, i64, char, char>, i64> PairVotes;

static void add_pair_votes(const vector<pair<i64, char>> &observed, i64 support, PairVotes &votes) {
    vector<pair<i64, char>> o = observed;
    sort(o.begin(), o.end());
    for (size_t i = 0; i < o.size(); i++)
        for (size_t j = i + 1; j < o.size(); j++) votes[{o[i].first, o[j].first, o[i].second, o[j].second}] += support;
}

// genotypes: (coordinate, (nt1, nt2)) in coordinate order; returns haplotypes and the phased coordinates
static void phase_from_votes(const vector<pair<i64, pair<char, char>>> &genotypes, const PairVotes &pair_votes, int min_support,
                             map<i64, pair<char, char>> &haplotypes, set<i64> &phased) {
    map<i64, pair<char, char>> g(genotypes.begin(), genotypes.end());
    vector<i64> het;
    for (auto &x : genotypes) if (x.second.first != x.second.second) het.push_back(x.first);
    sort(het.begin(), het.end());
    map<pair<i64, i64>, i64> votes;
    if (het.size() > 1) {
        set<i64> hs(het.begin(), het.end());
        for (auto &[k, support] : pair_votes) {
            auto [c1, c2, n1, n2] = k;
            if (!hs.count(c1) || !hs.count(c2)) continue;
            auto &g1 = g[c1], &g2 = g[c2];
            if ((n1 != g1.first && n1 != g1.second) || (n2 != g2.first && n2 != g2.second)) continue;
            bool same = (n1 == g1.first) == (n2 == g2.first);
            votes[{c1, c2}] += same ? support : -support;
        }
    }
    ParityUnionFind forest(het);
    vector<pair<pair<i64, i64>, i64>> sv(votes.begin(), votes.end());
    stable_sort(sv.begin(), sv.end(), [](const auto &x, const auto &y) { return llabs(x.second) > llabs(y.second); });
    for (auto &[cc, vote] : sv) if (llabs(vote) >= min_support) forest.unite(cc.first, cc.second, vote > 0 ? 0 : 1);
    vector<i64> roots; map<i64, vector<i64>> blocks;
    for (i64 c : het) {
        i64 r = forest.find(c).first;
        if (!blocks.count(r)) roots.push_back(r);
        blocks[r].push_back(c);
    }
    set<i64> main_block;
    size_t best = 0;
    for (i64 r : roots) if (blocks[r].size() > best) { best = blocks[r].size(); main_block = set<i64>(blocks[r].begin(), blocks[r].end()); }
    haplotypes.clear(); phased.clear();
    for (auto &[c, ab] : genotypes) {
        char a = ab.first, b = ab.second;
        if (a != b && main_block.count(c) && forest.find(c).second == 1) swap(a, b);
        haplotypes[c] = {a, b};
        if (ab.first == ab.second) phased.insert(c);
    }
    phased.insert(main_block.begin(), main_block.end());
}

// ------------------------------------------------------------------ call: loci with several bubbles
struct LocusData {
    int locus;
    vector<i64> coordinates;
    map<i64, vector<char>> alleles;
    set<i64> site_specific;
    vector<int> multi_snp;
    int n_real = 0, n_synthetic = 0;
    i64 site_base = 0, sa_base = 0, pl_base = 0;
};

static vector<i64> bubble_depth_vector(const Materialised &M, const Store &st, i64 coordinate, const vector<Carrier> &carriers) {
    struct Member { int bubble, path; set<i64> own; };
    vector<vector<Member>> groups;
    for (auto &cr : carriers) {
        const Bubble &b = M.bubbles[cr.bubble];
        set<i64> own;
        for (auto &s : b.snps) own.insert(b.coordinate(s.pos));
        bool placed = false;
        for (auto &group : groups) {
            for (auto &other : group) {
                const Bubble &ob = M.bubbles[other.bubble];
                set<i64> un = own;
                un.insert(other.own.begin(), other.own.end());
                un.erase(coordinate);
                bool differ = false;
                for (i64 c : un) {
                    char n1 = b.nucleotide(cr.path, c), n2 = ob.nucleotide(other.path, c);
                    if (n1 != NONE && n2 != NONE && n1 != n2) { differ = true; break; }
                }
                if (!differ) { group.push_back({cr.bubble, cr.path, own}); placed = true; break; }
            }
            if (placed) break;
        }
        if (!placed) groups.push_back({{cr.bubble, cr.path, own}});
    }
    vector<i64> total(st.S, 0);
    for (auto &group : groups) {
        vector<i64> mx(st.S, 0);
        for (auto &m : group) {
            size_t r = (size_t) (2 * M.bubbles[m.bubble].index + m.path) * st.S;
            for (int s = 0; s < st.S; s++) mx[s] = max<i64>(mx[s], st.counts[r + s]);
        }
        for (int s = 0; s < st.S; s++) total[s] += mx[s];
    }
    for (auto &x : total) x = min<i64>(x, 65535);
    return total;
}

// ---- repeat and low complexity scores (see the python script: low_complexity, repeat_scores)
static const int DUST_WINDOW = 64, DUST_LEVEL = 20;
struct RepeatScores { vector<i64> degree; vector<char> giant; double median_depth = 0; int max_overlaps = 50; };
static RepeatScores SCORES;

static double low_complexity(const string &seq, i64 center) {
    i64 b = min<i64>((i64) seq.size(), max<i64>(0, center - DUST_WINDOW / 2) + DUST_WINDOW);
    i64 a = max<i64>(0, b - DUST_WINDOW);
    int counts[64] = {0};
    i64 n = 0;
    auto code = [](char c) -> int {
        c = (char) toupper(c);
        return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
    };
    for (i64 i = a; i + 2 < b; i++) {
        int x = code(seq[i]), y = code(seq[i + 1]), z = code(seq[i + 2]);
        if (x < 0 || y < 0 || z < 0) continue;
        counts[x * 16 + y * 4 + z]++;
        n++;
    }
    if (n <= 1) return 0.0;
    i64 sum = 0;
    for (int c : counts) sum += (i64) c * (c - 1) / 2;
    return min(1.0, (double) sum / (double) (n - 1) / DUST_LEVEL);
}
static void repeat_scores(i64 degree, i64 depth, double &ovs, double &dps) {
    int reference = SCORES.max_overlaps > 0 ? SCORES.max_overlaps : 50;
    ovs = min(1.0, log1p((double) degree) / log1p((double) reference));
    dps = SCORES.median_depth <= 0 ? 0.0 : 1.0 - exp(-max(0.0, (double) depth / SCORES.median_depth - 1.0));
}
// RPT = max(OVS, DPS), 1 for the bubbles of a giant component (GC flag); LONG: locus longer than --max_locus_length
static string repeat_info(double ovs, double dps, double lc, bool giant = false, bool too_long = false) {
    char b[160];
    snprintf(b, sizeof(b), ";RPT=%.3f;OVS=%.3f;DPS=%.3f;LC=%.3f", giant ? 1.0 : max(ovs, dps), ovs, dps, lc);
    return string(b) + (giant ? ";GC" : "") + (too_long ? ";LONG" : "");
}
static bool too_long(i64 length, int max_locus_length) { return max_locus_length > 0 && length > max_locus_length; }
static double median_of(vector<i64> v) {
    if (v.empty()) return 0.0;
    sort(v.begin(), v.end());
    size_t n = v.size();
    return n % 2 ? (double) v[n / 2] : (double) (v[n / 2 - 1] + v[n / 2]) / 2;
}

static string info_field(double rank, bool has_rank, const array<int, 4> &meta, i64 cluster, i64 cluster_size, i64 n_real,
                         i64 n_synthetic, const vector<i64> &ids, i64 n_sites, i64 conflicts) {
    string m[4];
    for (int i = 0; i < 4; i++) m[i] = meta[i] < 0 ? "." : S(meta[i]);
    string r = (!has_rank || std::isnan(rank)) ? "." : py_float(rank);
    return "Ty=SNP;Rk=" + r + ";UL=" + m[0] + ";UR=" + m[1] + ";CL=" + m[2] + ";CR=" + m[3] + ";Genome=.;Sd=.;Cluster=" + S(cluster)
           + ";ClSize=" + S(cluster_size) + ";NB=" + S(n_real) + ";NX=" + S(n_synthetic) + ";BUB=" + join_ints(ids) + ";NSITES="
           + S(n_sites) + ";PC=" + S(conflicts);
}

static void write_alleles(Out &handle, const string &name, const string &header, const string &templ, const vector<i64> &coords,
                          const map<string, i64> &copies, const string &reference, const map<string, i64> *partial) {
    auto order = [&](const map<string, i64> &c) {
        vector<pair<string, i64>> v(c.begin(), c.end());
        stable_sort(v.begin(), v.end(), [&](const auto &x, const auto &y) {
            return make_tuple(-x.second, x.first != reference, x.first) < make_tuple(-y.second, y.first != reference, y.first);
        });
        return v;
    };
    vector<tuple<string, i64, string>> haps;
    for (auto &x : order(copies)) haps.push_back({x.first, x.second, "resolved"});
    if (partial) for (auto &x : order(*partial)) haps.push_back({x.first, x.second, "partial"});
    if (haps.empty()) haps.push_back({reference, 0, "reference"});
    for (size_t r = 0; r < haps.size(); r++) {
        auto &[hap, n, status] = haps[r];
        string seq = templ;
        for (size_t i = 0; i < coords.size() && i < hap.size(); i++) seq[coords[i]] = hap[i];
        string path = r == 0 ? name + "_higher_path" : name + "_lower_path_" + S(r);
        handle << ">" << path << "|haplotype_" << hap << "|copies_" << n << "|status_" << status << "|" << header << "\n" << seq << "\n";
    }
}

static string locus_sequence(const Materialised &M, const Locus &locus, const map<i64, char> &site_nucleotides) {
    vector<array<int, 5>> counts(locus.length, {0, 0, 0, 0, 0});
    vector<char> upper(locus.length, 0);
    for (int s : locus.bubbles) {
        const Bubble &b = M.bubbles[s];
        if (b.synthetic || b.full.empty()) continue;
        size_t L = b.full.size();
        auto ex = b.extent();
        i64 start = min(b.coordinate(ex.first), b.coordinate(ex.second));
        for (size_t j = 0; j < L; j++) {
            char ch = b.orient == 1 ? b.full[j] : b.full[L - 1 - j];
            uint8_t code = CODE[(uint8_t) ch];
            if (b.orient == -1) code = COMP[code];
            counts[start + j][min<int>(code, 4)]++;
            if ((uint8_t) ch < 'a') upper[start + j] = 1;
        }
    }
    string seq(locus.length, 'N');
    for (i64 i = 0; i < locus.length; i++) {
        auto &c = counts[i];
        if (c[0] + c[1] + c[2] + c[3] > 0) {
            int best = 0;
            for (int k = 1; k < 4; k++) if (c[k] > c[best]) best = k;
            seq[i] = DECODE[best];
        }
    }
    for (auto &[c, n] : site_nucleotides) seq[c] = n;
    for (i64 i = 0; i < locus.length; i++) if (!upper[i]) seq[i] = (char) tolower(seq[i]);
    return seq;
}

static const char *VCF_HEADER =
    "##fileformat=VCFv4.2\n"
    "##source=disco_haplotypes.py\n"
    "##INFO=<ID=Ty,Number=1,Type=String,Description=\"SNP, INS, DEL or .\">\n"
    "##INFO=<ID=Rk,Number=1,Type=Float,Description=\"SNP rank (best rank among the kissnp2 bubbles describing this site)\">\n"
    "##INFO=<ID=UL,Number=1,Type=Integer,Description=\"length of the unitig left (of the best ranked bubble)\">\n"
    "##INFO=<ID=UR,Number=1,Type=Integer,Description=\"length of the unitig right (of the best ranked bubble)\">\n"
    "##INFO=<ID=CL,Number=1,Type=Integer,Description=\"length of the contig left (of the best ranked bubble)\">\n"
    "##INFO=<ID=CR,Number=1,Type=Integer,Description=\"length of the contig right (of the best ranked bubble)\">\n"
    "##INFO=<ID=Genome,Number=1,Type=String,Description=\"Allele of the reference;for indel reference is . \">\n"
    "##INFO=<ID=Sd,Number=1,Type=Integer,Description=\"Reverse (-1) or Forward (1) Alignement\">\n"
    "##INFO=<ID=Cluster,Number=1,Type=Integer,Description=\"Locus (cluster) number, as in CHROM\">\n"
    "##INFO=<ID=ClSize,Number=1,Type=Integer,Description=\"Cluster size: number of kissnp2 bubble paths in the locus (2 per bubble)\">\n"
    "##INFO=<ID=NB,Number=1,Type=Integer,Description=\"Number of kissnp2 bubbles describing this site\">\n"
    "##INFO=<ID=NX,Number=1,Type=Integer,Description=\"Number of synthetic context bubbles describing this site\">\n"
    "##INFO=<ID=BUB,Number=.,Type=String,Description=\"Ids of the kissnp2 bubbles describing this site\">\n"
    "##INFO=<ID=NSITES,Number=1,Type=Integer,Description=\"Number of sites of the locus\">\n"
    "##INFO=<ID=PC,Number=1,Type=Integer,Description=\"Bubble placement conflicts in this locus (repeats, paralogs)\">\n"
    "##INFO=<ID=RPT,Number=1,Type=Float,Description=\"Repeat score, 0-1: max(OVS, DPS); 1 = repeated (multi-copy) region\">\n"
    "##INFO=<ID=OVS,Number=1,Type=Float,Description=\"Overlap score, 0-1: log(1 + bubbles overlapped by the bubbles of the site) / log(1 + --max_overlaps), capped at 1\">\n"
    "##INFO=<ID=DPS,Number=1,Type=Float,Description=\"Depth score, 0-1: 1 - exp(-(site depth / median bubble depth - 1)), 0 at or below the median\">\n"
    "##INFO=<ID=LC,Number=1,Type=Float,Description=\"Low complexity score, 0-1: DUST score of the 64 nt around the site / 20 (dustmasker level), capped at 1\">\n"
    "##INFO=<ID=GC,Number=0,Type=Flag,Description=\"Bubble of a component larger than --max_locus_bubbles (bubbles chained by repeats, genotyped one at a time): RPT set to 1\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype; | = phased with the other sites of the locus sharing its PS, / = not phased\">\n"
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set (position of its first site)\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sum of the allele depths\">\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths (lower bounds, see the documentation)\">\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";

// observations of one fragment in one locus, in order of first appearance; NONE = contradictory
typedef vector<pair<i64, char>> Observed;
static void observe(Observed &obs, i64 c, char n, i64 &conflicts) {
    for (auto &x : obs) if (x.first == c) {
        if (x.second != n) { x.second = NONE; conflicts++; }
        return;
    }
    obs.push_back({c, n});
}

static i64 call_multi_loci(vector<Locus> &loci, Materialised &M, const Store &st, const vector<string> &fact_files,
                           const vector<string> &site_files, const Args &args, Out &vcf, Out &hap_file, Out &locus_file, Out &fasta,
                           Out &alleles_fasta) {
    int S_ = st.S;
    vector<LocusData> data(loci.size());
    unordered_map<i64, int> by_id;                                    // bubble id -> slot (loci kept)
    vector<int> locus_data_of(loci.size(), -1);
    for (size_t li = 0; li < loci.size(); li++) {
        LocusData &d = data[li];
        Locus &locus = loci[li];
        d.locus = (int) li;
        for (auto &[c, carriers] : locus.sites) {
            d.coordinates.push_back(c);
            int first = -1;
            for (auto &x : carriers.v) for (auto &cr : x.second)
                if (first < 0 || M.bubbles[cr.bubble].index < M.bubbles[first].index) first = cr.bubble;
            char reference = M.bubbles[first].nucleotide(0, c);
            vector<char> others;
            for (auto &x : carriers.v) if (x.first != reference) others.push_back(x.first);
            sort(others.begin(), others.end());
            vector<char> al = {reference};
            al.insert(al.end(), others.begin(), others.end());
            d.alleles[c] = al;
            bool specific = false;
            for (auto &x : carriers.v) for (auto &cr : x.second) if (M.bubbles[cr.bubble].snps.size() == 1) specific = true;
            if (specific) d.site_specific.insert(c);
        }
        for (int s : locus.bubbles) {
            if (M.bubbles[s].snps.size() > 1) d.multi_snp.push_back(s);
            if (!M.bubbles[s].synthetic) d.n_real++;
            by_id[M.bubbles[s].id] = s;
        }
        d.n_synthetic = (int) locus.bubbles.size() - d.n_real;
    }
    // ---- global slots
    i64 site_base = 0, sa_base = 0, pl_base = 0;
    for (auto &d : data) {
        d.site_base = site_base; d.sa_base = sa_base; d.pl_base = pl_base;
        for (i64 c : d.coordinates) {
            i64 A = d.alleles[c].size();
            site_base++; sa_base += A; pl_base += A * (A + 1) / 2;
        }
    }
    i64 n_sites = site_base, n_sa = sa_base, n_pl = pl_base;
    vector<uint16_t> AD((size_t) n_sa * S_, 0), PL((size_t) n_pl * S_, 0);
    vector<int8_t> GT((size_t) n_sites * S_ * 2, -1);
    vector<char> PH((size_t) n_sites * S_, 0);
    vector<int32_t> PS((size_t) n_sites * S_, 0);
    vector<uint8_t> GQ((size_t) n_sites * S_, 0);
    vector<int> site_A(n_sites), site_locus(n_sites);
    vector<i64> site_sa(n_sites), site_pl(n_sites), site_coord(n_sites);
    for (size_t li = 0; li < data.size(); li++) {
        LocusData &d = data[li];
        i64 site = d.site_base, sa = d.sa_base, pl = d.pl_base;
        for (i64 c : d.coordinates) {
            auto &al = d.alleles[c];
            site_A[site] = (int) al.size(); site_sa[site] = sa; site_pl[site] = pl; site_coord[site] = c; site_locus[site] = (int) li;
            const SiteCarriers &carriers = loci[li].sites.at(c);
            for (char nt : al) {
                const vector<Carrier> &all = *carriers.get(nt);
                vector<Carrier> counted;
                for (auto &cr : all) if (M.bubbles[cr.bubble].snps.size() == 1) counted.push_back(cr);
                if (counted.empty()) counted = all;
                vector<i64> dv = bubble_depth_vector(M, st, c, counted);
                for (int s = 0; s < S_; s++) AD[(size_t) sa * S_ + s] = (uint16_t) dv[s];
                sa++;
            }
            site++;
            pl += (i64) al.size() * (al.size() + 1) / 2;
        }
    }
    // row of AD of (locus, coordinate, nucleotide), -1
    auto sa_row = [&](int li, i64 c, char n) -> i64 {
        const LocusData &d = data[li];
        auto it = lower_bound(d.coordinates.begin(), d.coordinates.end(), c);
        if (it == d.coordinates.end() || *it != c) return -1;
        i64 site = d.site_base + (it - d.coordinates.begin());
        const auto &al = d.alleles.at(c);
        for (size_t k = 0; k < al.size(); k++) if (al[k] == n) return site_sa[site] + (i64) k;
        return -1;
    };
    map<int, pair<string, bool>> fact_of_sample;
    for (auto &f : fact_files) fact_of_sample[sample_of_fact_file(f)] = {f, false};
    for (auto &f : site_files) fact_of_sample[sample_of_fact_file(f)] = {f, true};
    map<int, string> names;
    for (auto &[s, fb] : fact_of_sample) names[s] = read_set_name(fb.first);
    i64 total_conflicts = 0;

    // ---- sample by sample (independent: in parallel)
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : total_conflicts)
    for (int sample = 0; sample < S_; sample++) {
        vector<PairVotes> votes(loci.size());
        vector<char> has_votes(loci.size(), 0);
        auto fos = fact_of_sample.find(sample);
        if (fos != fact_of_sample.end()) {
            const string &file = fos->second.first;
            bool site_facts = fos->second.second;
            unordered_map<i64, i64> fact_depth;
            LineReader in(file);
            string line, text;
            i64 support;
            vector<pair<int, Observed>> per_locus;       // locus -> observations, in order of first appearance
            vector<vector<Token>> parts;
            vector<string> toks;
            while (in.next(line)) {
                line += '\n';
                if (!fact_line(line, text, support)) continue;
                per_locus.clear();
                auto locus_obs = [&](int li) -> Observed & {
                    for (auto &x : per_locus) if (x.first == li) return x.second;
                    per_locus.push_back({li, {}});
                    return per_locus.back().second;
                };
                if (site_facts) {
                    // '12h:101;15l:1;'
                    size_t b = 0, e = text.size();
                    while (b < e && isspace((unsigned char) text[b])) b++;
                    while (e > b && isspace((unsigned char) text[e - 1])) e--;
                    size_t p = b;
                    while (p <= e) {
                        size_t q = text.find(';', p);
                        if (q == string::npos || q > e) q = e;
                        if (q > p) {
                            // (\d+)([hl]):([01]+)$
                            size_t x = p; i64 id = 0;
                            while (x < q && isdigit((unsigned char) text[x])) id = id * 10 + (text[x++] - '0');
                            bool ok = x > p && x < q && (text[x] == 'h' || text[x] == 'l') && x + 1 < q && text[x + 1] == ':';
                            size_t ms = x + 2;
                            if (ok) {
                                if (ms >= q) ok = false;
                                for (size_t y = ms; ok && y < q; y++) if (text[y] != '0' && text[y] != '1') ok = false;
                            }
                            if (ok) {
                                auto it = by_id.find(id);
                                if (it != by_id.end()) {
                                    const Bubble &bb = M.bubbles[it->second];
                                    if (q - ms == bb.snps.size()) {
                                        int pi = text[x] == 'h' ? 0 : 1;
                                        Observed &obs = locus_obs(bb.locus);
                                        for (size_t k = 0; k < bb.snps.size(); k++) {
                                            if (text[ms + k] != '1') continue;
                                            i64 c = bb.coordinate(bb.snps[k].pos);
                                            observe(obs, c, bb.nucleotide(pi, c), total_conflicts);
                                        }
                                    }
                                }
                            }
                        }
                        p = q + 1;
                    }
                } else {
                    fact_parts(text, parts);
                    for (auto &part : parts) for (auto &t : part) {
                        auto it = by_id.find(t.id);
                        if (it == by_id.end()) continue;
                        const Bubble &bb = M.bubbles[it->second];
                        Observed &obs = locus_obs(bb.locus);
                        const set<i64> &specific = data[bb.locus].site_specific;
                        for (auto &snp : bb.snps) {
                            i64 c = bb.coordinate(snp.pos);
                            if (bb.snps.size() > 1 && specific.count(c)) continue;
                            observe(obs, c, bb.nucleotide(t.path, c), total_conflicts);
                        }
                    }
                }
                for (auto &[li, obs] : per_locus) {
                    Observed kept;
                    for (auto &x : obs) if (x.second != NONE) kept.push_back(x);
                    if (kept.empty()) continue;
                    for (auto &[c, n] : kept) {
                        i64 row = sa_row(li, c, n);
                        if (row >= 0) fact_depth[row] += support;
                    }
                    has_votes[li] = 1;
                    add_pair_votes(kept, support, votes[li]);
                }
            }
            for (auto &[row, value] : fact_depth) {
                uint16_t &a = AD[(size_t) row * S_ + sample];
                a = (uint16_t) max<i64>(a, min<i64>(value, 65535));
            }
        }
        // genotypes of every site of this sample
        vector<uint16_t> depths, pl;
        for (i64 site = 0; site < n_sites; site++) {
            int A = site_A[site];
            depths.resize(A);
            for (int a = 0; a < A; a++) depths[a] = AD[(size_t) (site_sa[site] + a) * S_ + sample];
            int gt[2], gq;
            call_genotype(depths.data(), A, args.min_depth, gt, pl, gq);
            GT[((size_t) site * S_ + sample) * 2] = (int8_t) gt[0];
            GT[((size_t) site * S_ + sample) * 2 + 1] = (int8_t) gt[1];
            GQ[(size_t) site * S_ + sample] = (uint8_t) gq;
            for (size_t g = 0; g < pl.size(); g++) PL[(size_t) (site_pl[site] + g) * S_ + sample] = pl[g];
        }
        // phasing
        for (size_t li = 0; li < data.size(); li++) {
            const LocusData &d = data[li];
            size_t ns = d.coordinates.size();
            int n_het = 0;
            i64 first_called = -1;
            for (size_t i = 0; i < ns; i++) {
                size_t k = ((size_t) (d.site_base + i) * S_ + sample) * 2;
                if (GT[k] >= 0) { if (first_called < 0) first_called = d.coordinates[i]; if (GT[k] != GT[k + 1]) n_het++; }
            }
            if (n_het < 2) {
                if (ns > 1)
                    for (size_t i = 0; i < ns; i++) {
                        size_t site = d.site_base + i;
                        if (GT[((size_t) site * S_ + sample) * 2] >= 0) {
                            PH[(size_t) site * S_ + sample] = 1;
                            PS[(size_t) site * S_ + sample] = (int32_t) (first_called + 1);
                        }
                    }
                continue;
            }
            vector<pair<i64, pair<char, char>>> called;
            for (size_t i = 0; i < ns; i++) {
                size_t k = ((size_t) (d.site_base + i) * S_ + sample) * 2;
                if (GT[k] >= 0) {
                    const auto &al = d.alleles.at(d.coordinates[i]);
                    called.push_back({d.coordinates[i], {al[GT[k]], al[GT[k + 1]]}});
                }
            }
            if (called.empty()) continue;
            PairVotes local;
            const PairVotes *pv = &votes[li];
            if (!has_votes[li]) {
                for (int s : d.multi_snp) {
                    const Bubble &bb = M.bubbles[s];
                    for (int pi = 0; pi < 2; pi++) {
                        Observed obs;
                        for (auto &snp : bb.snps) {
                            i64 c = bb.coordinate(snp.pos);
                            char n = bb.nucleotide(pi, c);
                            bool found = false;
                            for (auto &x : obs) if (x.first == c) { x.second = n; found = true; }
                            if (!found) obs.push_back({c, n});
                        }
                        add_pair_votes(obs, st.counts[(size_t) (2 * bb.index + pi) * S_ + sample], local);
                    }
                }
                pv = &local;
            }
            map<i64, pair<char, char>> haplotypes;
            set<i64> phased;
            phase_from_votes(called, *pv, args.min_support, haplotypes, phased);
            i64 phase_set = phased.empty() ? 0 : *phased.begin() + 1;
            for (size_t i = 0; i < ns; i++) {
                i64 c = d.coordinates[i];
                auto it = haplotypes.find(c);
                if (it == haplotypes.end()) continue;
                const auto &al = d.alleles.at(c);
                size_t site = d.site_base + i;
                int a = (int) (find(al.begin(), al.end(), it->second.first) - al.begin());
                int b = (int) (find(al.begin(), al.end(), it->second.second) - al.begin());
                GT[((size_t) site * S_ + sample) * 2] = (int8_t) a;
                GT[((size_t) site * S_ + sample) * 2 + 1] = (int8_t) b;
                if (phased.count(c)) { PH[(size_t) site * S_ + sample] = 1; PS[(size_t) site * S_ + sample] = (int32_t) phase_set; }
            }
        }
    }

    // ---- output
    i64 n_multiallelic = 0, sum_conflicts = 0;
    for (size_t li = 0; li < data.size(); li++) {
        LocusData &d = data[li];
        Locus &locus = loci[li];
        i64 number = (i64) li + 1;
        locus.name = "locus_" + S(number);
        sum_conflicts += locus.conflicts;
        string consensus = locus_sequence(M, locus, {});
        for (size_t i = 0; i < d.coordinates.size(); i++) {
            i64 c = d.coordinates[i];
            i64 site = d.site_base + i;
            auto &al = d.alleles[c];
            int A = (int) al.size();
            n_multiallelic += A > 2;
            vector<int> involved;
            for (auto &x : locus.sites.at(c).v) for (auto &cr : x.second)
                if (find(involved.begin(), involved.end(), cr.bubble) == involved.end()) involved.push_back(cr.bubble);
            sort(involved.begin(), involved.end(), [&](int x, int y) { return M.bubbles[x].index < M.bubbles[y].index; });
            vector<int> real;
            for (int s : involved) if (!M.bubbles[s].synthetic) real.push_back(s);
            if (real.empty()) real = involved;
            int best = real[0];
            auto rank_of = [&](int s) { return M.bubbles[s].has_rank() ? M.bubbles[s].rank : -1.0; };
            for (int s : real) if (rank_of(s) > rank_of(best)) best = s;
            vector<i64> ids;
            i64 n_r = 0, n_x = 0;
            for (int s : involved) { if (M.bubbles[s].synthetic) n_x++; else { n_r++; ids.push_back(M.bubbles[s].id); } }
            sort(ids.begin(), ids.end());
            string info = info_field(M.bubbles[best].rank, M.bubbles[best].has_rank(), M.bubbles[best].meta, number, 2 * d.n_real, n_r, n_x,
                                     ids, d.coordinates.size(), locus.conflicts);
            i64 degree = 0, depth = 0;
            for (int s : involved) if (!M.bubbles[s].synthetic) degree = max(degree, SCORES.degree[M.bubbles[s].index]);
            for (int a = 0; a < A; a++) for (int s = 0; s < S_; s++) depth += AD[(size_t) (site_sa[site] + a) * S_ + s];
            double ovs, dps;
            repeat_scores(degree, depth, ovs, dps);
            info += repeat_info(ovs, dps, low_complexity(consensus, c), false, too_long(locus.length, args.max_locus_length));
            string alt;
            for (int k = 1; k < A; k++) { if (k > 1) alt += ","; alt += al[k]; }
            vcf << locus.name << "\t" << (c + 1) << "\t" << locus.name << "_" << (c + 1) << "\t" << string(1, al[0]) << "\t" << alt
                << "\t.\t.\t" << info << "\tGT:PS:DP:AD:GQ:PL";
            for (int s = 0; s < S_; s++) {
                vector<i64> ad(A), pl(A * (A + 1) / 2);
                for (int a = 0; a < A; a++) ad[a] = AD[(size_t) (site_sa[site] + a) * S_ + s];
                for (size_t g = 0; g < pl.size(); g++) pl[g] = PL[(size_t) (site_pl[site] + g) * S_ + s];
                int gt[2] = {GT[((size_t) site * S_ + s) * 2], GT[((size_t) site * S_ + s) * 2 + 1]};
                vcf << "\t" << format_field(gt, PH[(size_t) site * S_ + s], PS[(size_t) site * S_ + s], ad, GQ[(size_t) site * S_ + s], pl);
            }
            vcf << "\n";
        }
        map<string, i64> seen, partial;
        string positions;
        for (size_t i = 0; i < d.coordinates.size(); i++) { if (i) positions += ","; positions += S(d.coordinates[i] + 1); }
        for (int s = 0; s < S_; s++) {
            string h1, h2, a1, a2;
            i64 n_missing = 0, n_unphased = 0;
            for (size_t i = 0; i < d.coordinates.size(); i++) {
                i64 c = d.coordinates[i];
                size_t site = d.site_base + i;
                int a = GT[((size_t) site * S_ + s) * 2], b = GT[((size_t) site * S_ + s) * 2 + 1];
                auto &al = d.alleles[c];
                if (a < 0) { h1 += 'N'; h2 += 'N'; a1 += 'N'; a2 += 'N'; n_missing++; }
                else if (!PH[site * S_ + s] && a != b && d.coordinates.size() > 1) {
                    h1 += '?'; h2 += '?'; n_unphased++;
                    char code = iupac({al[a], al[b]});
                    a1 += code; a2 += code;
                } else { h1 += al[a]; h2 += al[b]; a1 += al[a]; a2 += al[b]; }
            }
            string status;
            if (n_missing == (i64) d.coordinates.size()) status = "missing";
            else if (n_missing || n_unphased) {
                status = "partial";
                partial[min(a1, a2)]++; partial[max(a1, a2)]++;
            } else {
                status = "resolved";
                if (h2 < h1) swap(h1, h2);
                seen[h1]++; seen[h2]++;
            }
            auto nm = names.find(s);
            hap_file << locus.name << "\t" << d.coordinates.size() << "\t" << positions << "\tG" << (s + 1) << "\t"
                     << (nm == names.end() ? string(".") : nm->second) << "\t" << status << "\t" << h1 << "\t" << h2 << "\t"
                     << n_missing << "\t" << n_unphased << "\n";
        }
        size_t max_alleles = 0;
        for (auto &[c, al] : d.alleles) max_alleles = max(max_alleles, al.size());
        string haps;
        for (auto &[h, n] : seen) { if (!haps.empty()) haps += ","; haps += h + ":" + S(n); }
        locus_file << locus.name << "\t" << locus.length << "\t" << d.n_real << "\t" << d.n_synthetic << "\t" << d.coordinates.size()
                   << "\t" << max_alleles << "\t" << seen.size() << "\t" << (haps.empty() ? string(".") : haps) << "\t" << locus.conflicts
                   << "\n";
        string header = "length_" + S(locus.length) + "|n_sites_" + S(d.coordinates.size()) + "|positions_" + positions;
        map<i64, char> iup;
        for (auto &[c, al] : d.alleles) iup[c] = iupac(al);
        fasta << ">" << locus.name << "|" << header << "\n" << locus_sequence(M, locus, iup) << "\n";
        string reference;
        for (i64 c : d.coordinates) reference += d.alleles[c][0];
        write_alleles(alleles_fasta, locus.name, header, consensus, d.coordinates, seen, reference, &partial);
    }
    log_msg("[call] " + S(data.size()) + " loci with several bubbles: " + S(n_sites) + " sites, " + S(n_multiallelic)
            + " with more than two alleles, " + S(sum_conflicts) + " placement conflicts, " + S(total_conflicts)
            + " contradictory observations in facts");
    return (i64) data.size();
}

// ------------------------------------------------------------------ call: isolated bubbles
struct SingleFasta { string name, pos_text; vector<i64> coords; string hap_h; i64 n_h; string hap_l; i64 n_l; };

static unordered_map<i64, SingleFasta> call_single_bubbles(const Store &st, const vector<int> &indices, i64 first_number, const Args &args,
                                                           Out &vcf, Out &hap_file, const map<int, string> &names, Out &locus_file) {
    int S_ = st.S;
    i64 number = first_number, n_sites = 0;
    unordered_map<i64, SingleFasta> for_fasta;
    vector<uint16_t> pl;
    for (int index : indices) {
        vector<int> positions = st.snp_positions(index);
        if (positions.empty()) continue;
        const string &higher = st.paths[2 * index], &lower = st.paths[2 * index + 1];
        string name = "locus_" + S(number);
        number++;
        i64 offset = st.left[2 * index];
        bool phased = positions.size() > 1;
        i64 phase_set = phased ? offset + positions[0] + 1 : 0;
        vector<array<int, 2>> gts(S_);
        string sample_columns;
        for (int s = 0; s < S_; s++) {
            uint16_t dp[2] = {st.counts[(size_t) (2 * index) * S_ + s], st.counts[(size_t) (2 * index + 1) * S_ + s]};
            int gt[2], gq;
            call_genotype(dp, 2, args.min_depth, gt, pl, gq);
            gts[s] = {gt[0], gt[1]};
            if (s) sample_columns += "\t";
            sample_columns += format_field(gt, phased, phase_set, {dp[0], dp[1]}, gq, {pl[0], pl[1], pl[2]});
        }
        double rank = st.ranks[index];
        string info = info_field(rank, !std::isnan(rank), st.meta[index], number - 1, 2, 1, 0, {st.ids[index]}, positions.size(), 0);
        i64 depth = 0;
        for (int s = 0; s < S_; s++) depth += (i64) st.counts[(size_t) (2 * index) * S_ + s] + st.counts[(size_t) (2 * index + 1) * S_ + s];
        double ovs, dps;
        repeat_scores(SCORES.degree[index], depth, ovs, dps);
        string extended(st.flens[2 * index], 'N');
        for (int k = 0; k < st.flens[2 * index]; k++) extended[k] = DECODE[st.frow(2 * index)[k]];
        bool giant = SCORES.giant[index], long_locus = too_long(st.left[2 * index] + st.lens(2 * index) + st.right[2 * index], args.max_locus_length);
        for (int p : positions) {
            i64 pos = offset + p + 1;
            double lc = low_complexity(extended, (i64) st.fstart[2 * index] + p);
            vcf << name << "\t" << pos << "\t" << name << "_" << pos << "\t" << string(1, higher[p]) << "\t" << string(1, lower[p])
                << "\t.\t.\t" << info << repeat_info(ovs, dps, lc, giant, long_locus) << "\tGT:PS:DP:AD:GQ:PL\t" << sample_columns << "\n";
        }
        n_sites += positions.size();
        string hap_h, hap_l;
        for (int p : positions) { hap_h += higher[p]; hap_l += lower[p]; }
        i64 n_h = 0, n_l = 0;
        for (auto &g : gts) if (g[0] >= 0) { n_h += (g[0] == 0) + (g[1] == 0); n_l += (g[0] == 1) + (g[1] == 1); }
        vector<string> haplotypes;
        if (n_h) haplotypes.push_back(hap_h + ":" + S(n_h));
        if (n_l) haplotypes.push_back(hap_l + ":" + S(n_l));
        string pos_text;
        for (size_t i = 0; i < positions.size(); i++) { if (i) pos_text += ","; pos_text += S(offset + positions[i] + 1); }
        if (positions.size() > 1) {
            for (int s = 0; s < S_; s++) {
                auto nm = names.find(s);
                string rs = nm == names.end() ? "." : nm->second;
                hap_file << name << "\t" << positions.size() << "\t" << pos_text << "\tG" << (s + 1) << "\t" << rs << "\t";
                if (gts[s][0] < 0) {
                    string nn(positions.size(), 'N');
                    hap_file << "missing\t" << nn << "\t" << nn << "\t" << positions.size() << "\t0\n";
                } else {
                    hap_file << "resolved\t" << (gts[s][0] == 0 ? hap_h : hap_l) << "\t" << (gts[s][1] == 0 ? hap_h : hap_l) << "\t0\t0\n";
                }
            }
        }
        sort(haplotypes.begin(), haplotypes.end());
        string hs;
        for (size_t i = 0; i < haplotypes.size(); i++) { if (i) hs += ","; hs += haplotypes[i]; }
        locus_file << name << "\t" << (st.left[2 * index] + st.lens(2 * index) + st.right[2 * index]) << "\t1\t0\t" << positions.size()
                   << "\t2\t" << haplotypes.size() << "\t" << (hs.empty() ? string(".") : hs) << "\t0\n";
        SingleFasta sf{name, pos_text, {}, hap_h, n_h, hap_l, n_l};
        for (int p : positions) sf.coords.push_back(offset + p);
        for_fasta[st.ids[index]] = sf;
    }
    log_msg("[call] " + S(number - first_number) + " isolated bubbles: " + S(n_sites) + " sites");
    return for_fasta;
}

static void write_single_fastas(const string &file, const unordered_map<i64, SingleFasta> &for_fasta, Out &fasta, Out &alleles_fasta) {
    LineReader in(file);
    string line, header;
    bool has_header = false;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '>') { header = line; has_header = true; continue; }
        if (has_header && header.compare(0, 17, ">SNP_higher_path_") == 0) {
            bool snp, higher; i64 id;
            if (header_match(header, snp, higher, id)) {
                auto it = for_fasta.find(id);
                if (it != for_fasta.end()) {
                    const SingleFasta &r = it->second;
                    string templ = line;
                    while (!templ.empty() && (templ.back() == '\r' || templ.back() == '\n')) templ.pop_back();
                    string lh = "length_" + S(templ.size()) + "|n_sites_" + S(r.coords.size()) + "|positions_" + r.pos_text;
                    string seq = templ;
                    for (size_t i = 0; i < r.coords.size(); i++) seq[r.coords[i]] = iupac({r.hap_h[i], r.hap_l[i]});
                    fasta << ">" << r.name << "|" << lh << "\n" << seq << "\n";
                    map<string, i64> copies;
                    if (r.n_h) copies[r.hap_h] = r.n_h;
                    if (r.n_l) copies[r.hap_l] = r.n_l;
                    write_alleles(alleles_fasta, r.name, lh, templ, r.coords, copies, r.hap_h, nullptr);
                }
            }
        }
        has_header = false;
    }
}

static void load_full_paths(const string &file, Materialised &M) {
    LineReader in(file);
    string line, header;
    bool has_header = false;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '>') { header = line; has_header = true; continue; }
        if (has_header && header.compare(0, 17, ">SNP_higher_path_") == 0) {
            bool snp, higher; i64 id;
            if (header_match(header, snp, higher, id)) {
                auto it = M.slot_of_id.find(id);
                if (it != M.slot_of_id.end() && !M.bubbles[it->second].synthetic) {
                    string s = line;
                    while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back();
                    M.bubbles[it->second].full = s;
                }
            }
        }
        has_header = false;
    }
}

static void call(const Args &args) {
    SyntheticMap parents = read_synthetic_map(args.map);
    unordered_set<i64> synthetic_ids(parents.order.begin(), parents.order.end());
    Store st = parse_store(args.coherent, true, nullptr, args.max_flank);
    if (!args.uncoherent.empty() && !synthetic_ids.empty() && file_exists(args.uncoherent)) {
        Store extra = parse_store(args.uncoherent, true, &synthetic_ids, args.max_flank);
        if (extra.n) {
            vector<char> keep(extra.n);
            bool any = false;
            for (int i = 0; i < extra.n; i++) { keep[i] = st.index_of(extra.ids[i]) < 0; any |= keep[i]; }
            if (any) {
                if (extra.S != st.S) fatal("ERROR: the coherent and uncoherent files have different numbers of read sets");
                extend_store(st, extra, keep);
            }
        }
    }
    if (parents.order.empty())
        log_msg("[call] WARNING: no synthetic bubbles (no 'augment' step before kissreads2). Reads carrying\n"
                "[call]          other alleles at close SNPs are under-counted by kissreads2: genotypes and\n"
                "[call]          haplotypes of loci with several SNPs in one k-window are not reliable.");
    vector<string> fact_files, site_files;
    for (auto &f : args.phased) if (file_exists(f)) fact_files.push_back(f);
    for (auto &f : args.sites) if (file_exists(f)) site_files.push_back(f);
    if (site_files.empty())
        log_msg("[call] WARNING: no phased_sites files (kissreads2 -phasing_sites): reads mapping a multi-SNP bubble report\n"
                "[call]          only some of its sites, so some heterozygous sites may stay unphased.");
    if (fact_files.empty()) log_msg("[call] WARNING: no phased facts (kissreads2 -phasing): sites of different bubbles stay unphased.");

    Edges edges;
    for (i64 sid : parents.order) {
        auto &pp = parents.parent[sid];
        int s = st.index_of(sid), p = st.index_of(pp.first);
        if (s >= 0 && p >= 0) {
            if (p < s) { edges.b1.push_back(p); edges.b2.push_back(s); edges.shift.push_back(pp.second); }
            else { edges.b1.push_back(s); edges.b2.push_back(p); edges.shift.push_back(-pp.second); }
            edges.rel.push_back(1);
        }
    }
    View v{0, &st, {}};
    for (int i = 0; i < st.n; i++) if (!synthetic_ids.count(st.ids[i])) v.index.push_back(i);
    v.n = (int) v.index.size();
    Edges se = sequence_edges(v, args.seed_size, args.max_mismatches, args.min_overlap, args.max_divergence, args.seed_sampling);
    vector<i64> degree; vector<char> repeat_view;
    se = remove_repeats(se, v.n, args.max_overlaps, degree, repeat_view);
    for (size_t i = 0; i < se.size(); i++) {
        edges.b1.push_back(v.index[se.b1[i]]); edges.b2.push_back(v.index[se.b2[i]]);
        edges.shift.push_back(se.shift[i]); edges.rel.push_back(se.rel[i]);
    }
    Edges fe = fact_edges(fact_files, st);
    vector<char> kind(edges.size(), PARENT_EDGE);
    for (size_t i = edges.size() - se.size(); i < edges.size(); i++) kind[i] = SEQUENCE_EDGE;
    kind.resize(edges.size() + fe.size(), READ_EDGE);
    edges.b1.insert(edges.b1.end(), fe.b1.begin(), fe.b1.end()); edges.b2.insert(edges.b2.end(), fe.b2.begin(), fe.b2.end());
    edges.shift.insert(edges.shift.end(), fe.shift.begin(), fe.shift.end()); edges.rel.insert(edges.rel.end(), fe.rel.begin(), fe.rel.end());
    // the repeats are kept as isolated bubbles (no edge, their RPT is 1), their synthetic bubbles are removed
    vector<char> excluded(st.n, 0), repeat(st.n, 0);
    unordered_set<i64> repeat_ids;
    for (int i = 0; i < v.n; i++) if (repeat_view[i]) { repeat[v.index[i]] = 1; repeat_ids.insert(st.ids[v.index[i]]); }
    for (i64 sid : parents.order) {
        if (!repeat_ids.count(parents.parent[sid].first)) continue;
        int index = st.index_of(sid);
        if (index >= 0) excluded[index] = 1;
    }
    vector<i64> contradictions;
    {
        Edges kept;
        vector<char> kept_kind;
        for (size_t i = 0; i < edges.size(); i++) {
            if (excluded[edges.b1[i]] || excluded[edges.b2[i]] || repeat[edges.b1[i]] || repeat[edges.b2[i]]) continue;
            kept.b1.push_back(edges.b1[i]); kept.b2.push_back(edges.b2[i]); kept.shift.push_back(edges.shift[i]); kept.rel.push_back(edges.rel[i]);
            kept_kind.push_back(kind[i]);
        }
        // read links of multi-mapping reads, and with --max_locus_length the edges making a locus longer than
        // one ddRAD fragment or contradicting a placement, are removed
        EdgeRemovals removed;
        vector<char> keep = limit_locus_edges(kept, kept_kind, st, args.max_locus_length, args.min_overlap, args.max_mismatches,
                                              args.max_divergence, removed, contradictions);
        log_msg("[facts] " + S((i64) count(kept_kind.begin(), kept_kind.end(), READ_EDGE)) + " read links: " + S(removed.read_disagree)
                + " removed, the sequences of their bubbles disagree (multi-mapping reads)");
        if (args.max_locus_length > 0)
            log_msg("[loci] --max_locus_length " + S(args.max_locus_length) + ": removed " + S(removed.sequence_long)
                    + " sequence overlaps and " + S(removed.read_long) + " read links making a longer locus, "
                    + S(removed.sequence_contradiction) + " sequence overlaps and " + S(removed.read_contradiction)
                    + " read links contradicting a placement (counted as placement conflicts)");
        edges = Edges();
        for (size_t i = 0; i < kept.size(); i++) {
            if (!keep[i]) continue;
            edges.b1.push_back(kept.b1[i]); edges.b2.push_back(kept.b2[i]); edges.shift.push_back(kept.shift[i]); edges.rel.push_back(kept.rel[i]);
        }
        // repeat scores of the sites: overlapped bubbles, depth excess (median over the kissnp2 bubbles)
        SCORES.max_overlaps = args.max_overlaps;
        SCORES.degree.assign(st.n, 0);
        vector<i64> bubble_depth;
        for (int i = 0; i < v.n; i++) {
            int index = v.index[i];
            SCORES.degree[index] = degree[i];
            i64 d = 0;
            for (int s = 0; s < st.S; s++) d += (i64) st.counts[(size_t) (2 * index) * st.S + s] + st.counts[(size_t) (2 * index + 1) * st.S + s];
            bubble_depth.push_back(d);
        }
        SCORES.median_depth = median_of(bubble_depth);
        Out rep(args.out + "_repeats.tsv");
        rep << "bubble_id\toverlapped_bubbles\n";
        for (int i = 0; i < v.n; i++) if (repeat_view[i]) rep << st.ids[v.index[i]] << "\t" << degree[i] << "\n";
    }
    Materialised M;
    vector<char> is_synthetic(st.n);
    for (int i = 0; i < st.n; i++) is_synthetic[i] = synthetic_ids.count(st.ids[i]) > 0;
    SCORES.giant.assign(st.n, 0);
    vector<Locus> loci = multi_bubble_loci(st, edges, args.max_locus_bubbles, false, M, &is_synthetic, nullptr, &SCORES.giant);
    for (i64 index : contradictions)               // edges removed for contradicting the placement of their locus
        if (M.slot_of[index] >= 0 && M.bubbles[M.slot_of[index]].locus >= 0) loci[M.bubbles[M.slot_of[index]].locus].conflicts++;
    for (auto &b : M.bubbles) {
        auto it = parents.parent.find(b.id);
        if (it != parents.parent.end()) { b.parent = it->second.first; b.synthetic = true; }
    }
    load_full_paths(args.coherent, M);
    for (size_t i = 0; i < loci.size(); i++) loci[i].name = "locus_" + S(i + 1);
    vector<int> singles;
    i64 n_dropped = 0;
    for (int i = 0; i < st.n; i++) {
        bool in_locus = M.slot_of[i] >= 0, syn = synthetic_ids.count(st.ids[i]);
        if (excluded[i]) continue;
        if (!in_locus && !syn) singles.push_back(i);
        if (!in_locus && syn) n_dropped++;
    }
    if (n_dropped) log_msg("[call] " + S(n_dropped) + " synthetic bubbles without their parent bubble were dropped");
    if (args.min_sites > 1) {
        vector<Locus> kept;
        vector<int> new_index(loci.size(), -1);
        for (size_t i = 0; i < loci.size(); i++)
            if ((int) loci[i].sites.size() >= args.min_sites) { new_index[i] = (int) kept.size(); kept.push_back(loci[i]); }
        for (auto &b : M.bubbles) b.locus = b.locus >= 0 ? new_index[b.locus] : -1;
        loci.swap(kept);
        singles.clear();
    }

    Out vcf(args.out + ".vcf"), hap_file(args.out + ".tsv"), locus_file(args.out + "_loci.tsv"), fasta(args.out + "_loci.fa"),
        alleles_fasta(args.out + "_alleles.fa");
    {
        string header = VCF_HEADER;
        if (args.max_locus_length > 0) {
            string line = "##INFO=<ID=LONG,Number=0,Type=Flag,Description=\"Locus longer than --max_locus_length (" + S(args.max_locus_length) + " bp): longer than one ddRAD fragment sequenced by the reads (2 x read length - k), probably a chimera (paralogs, repeats)\">\n";
            header.insert(header.find("##FORMAT=<ID=GT,"), line);
        }
        vcf << header;
    }
    for (size_t i = 0; i < loci.size(); i++) vcf << "##contig=<ID=locus_" << (i64) (i + 1) << ",length=" << loci[i].length << ">\n";
    i64 number = (i64) loci.size() + 1;
    for (int index : singles)
        if (!st.snp_positions(index).empty()) {
            vcf << "##contig=<ID=locus_" << number << ",length=" << (st.left[2 * index] + st.lens(2 * index) + st.right[2 * index]) << ">\n";
            number++;
        }
    vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int s = 0; s < st.S; s++) vcf << "\tG" << (s + 1);
    vcf << "\n";
    hap_file << "locus\tn_sites\tpositions\tsample\tread_set\tstatus\thaplotype_1\thaplotype_2\tn_missing\tn_unphased\n";
    locus_file << "locus\tlength\tn_bubbles\tn_synthetic\tn_sites\tmax_alleles\tn_haplotypes\thaplotypes\tplacement_conflicts\n";
    i64 n_multi = call_multi_loci(loci, M, st, fact_files, site_files, args, vcf, hap_file, locus_file, fasta, alleles_fasta);
    if (!singles.empty()) {
        map<int, string> names;
        for (auto &f : fact_files) names[sample_of_fact_file(f)] = read_set_name(f);
        for (auto &f : site_files) names[sample_of_fact_file(f)] = read_set_name(f);
        auto for_fasta = call_single_bubbles(st, singles, n_multi + 1, args, vcf, hap_file, names, locus_file);
        write_single_fastas(args.coherent, for_fasta, fasta, alleles_fasta);
    }
}

// ------------------------------------------------------------------ strip
static void strip(const Args &args) {
    SyntheticMap m = read_synthetic_map(args.map);
    LineReader in(args.input);
    Out out(args.output);
    string line;
    bool keep = true;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '>') {
            keep = true;
            // >SNP_\w+_path_(\d+)
            if (line.compare(0, 5, ">SNP_") == 0) {
                size_t e = 5;
                while (e < line.size() && (isalnum((unsigned char) line[e]) || line[e] == '_')) e++;
                string word = line.substr(5, e - 5);
                // the last "_path_<digits>" inside the word run, preceded by at least one character
                size_t best = string::npos;
                for (size_t k = 1; k + 6 < word.size() + 1; k++)
                    if (word.compare(k, 6, "_path_") == 0 && k + 6 < word.size() && isdigit((unsigned char) word[k + 6])) best = k;
                if (best != string::npos) {
                    i64 id = 0;
                    for (size_t d = best + 6; d < word.size() && isdigit((unsigned char) word[d]); d++) id = id * 10 + (word[d] - '0');
                    keep = !m.parent.count(id);
                }
            }
        }
        if (keep) { out << line; if (in.ends_with_newline) out << '\n'; }
    }
}

// ------------------------------------------------------------------ main
static void usage() {
    cerr << "usage: disco_haplotypes {augment,call,strip} [options]\n"
            "  augment -i kissnp2.fa -o augmented.fa -m synthetic_map.tsv [--max_contexts 32] [--seed_window 31] [placement options]\n"
            "  call    -c coherent.fa [-u uncoherent.fa] [-m synthetic_map.tsv] [-p phased_alleles...] [-s phased_sites...] -o prefix\n"
            "          [--min_depth 3] [--min_support 2] [--min_sites 1] [--max_locus_length 0] [placement options]\n"
            "  strip   -i in.fa -o out.fa -m synthetic_map.tsv\n"
            "  placement options: --seed_size 16 --max_mismatches 1 --min_overlap 25 --max_flank 1000 --max_divergence 0.01\n"
            "                     --seed_sampling 8 --max_locus_bubbles 500 --max_overlaps 50 (0: no limit)\n"
            "                     --relaxed_placement (4 mismatches + 2% of the overlap, as the previous versions)\n"
            "  (same methods and outputs as scripts/disco_haplotypes.py)\n";
}

int main(int argc, char **argv) {
    init_tables();
    if (argc < 2) { usage(); return 2; }
    Args a;
    a.command = argv[1];
    if (a.command == "-h" || a.command == "--help") { usage(); return 0; }
    auto need = [&](int &i) -> string {
        if (i + 1 >= argc) { cerr << "ERROR: option " << argv[i] << " needs a value\n"; exit(2); }
        return argv[++i];
    };
    for (int i = 2; i < argc; i++) {
        string o = argv[i];
        auto many = [&](vector<string> &v) { while (i + 1 < argc && argv[i + 1][0] != '-') v.push_back(argv[++i]); };
        if (o == "-h" || o == "--help") { usage(); return 0; }
        else if (o == "-i" || o == "--input") a.input = need(i);
        else if (o == "-o" || o == "--output" || o == "--out") { if (a.command == "call") a.out = need(i); else a.output = need(i); }
        else if (o == "-m" || o == "--map") a.map = need(i);
        else if (o == "-c" || o == "--coherent") a.coherent = need(i);
        else if (o == "-u" || o == "--uncoherent") a.uncoherent = need(i);
        else if (o == "-p" || o == "--phased") many(a.phased);
        else if (o == "-s" || o == "--sites") many(a.sites);
        else if (o == "--max_contexts") a.max_contexts = stoi(need(i));
        else if (o == "--seed_window") a.seed_window = stoi(need(i));
        else if (o == "--min_depth") a.min_depth = stoi(need(i));
        else if (o == "--min_support") a.min_support = stoi(need(i));
        else if (o == "--min_sites") a.min_sites = stoi(need(i));
        else if (o == "--seed_size") a.seed_size = stoi(need(i));
        else if (o == "--max_mismatches") a.max_mismatches = stoi(need(i));
        else if (o == "--min_overlap") a.min_overlap = stoi(need(i));
        else if (o == "--max_flank") a.max_flank = stoi(need(i));
        else if (o == "--max_divergence") a.max_divergence = stod(need(i));
        else if (o == "--seed_sampling") a.seed_sampling = stoi(need(i));
        else if (o == "--max_locus_bubbles") a.max_locus_bubbles = stoi(need(i));
        else if (o == "--relaxed_placement") a.relaxed_placement = true;
        else if (o == "--max_overlaps") a.max_overlaps = stoi(need(i));
        else if (o == "--max_locus_length") a.max_locus_length = stoi(need(i));
        else { cerr << "ERROR: unknown option " << o << "\n"; usage(); return 2; }
    }
    if (a.max_flank > MAX_FLANK) { cerr << "ERROR: --max_flank cannot exceed " << MAX_FLANK << "\n"; return 2; }
    if (a.max_mismatches < 0) a.max_mismatches = a.relaxed_placement ? RELAXED_MISMATCHES : STRICT_MISMATCHES;
    if (a.max_divergence < 0) a.max_divergence = a.relaxed_placement ? RELAXED_DIVERGENCE : STRICT_DIVERGENCE;
    if (a.command == "augment") {
        if (a.input.empty() || a.output.empty() || a.map.empty()) { usage(); return 2; }
        augment(a);
    } else if (a.command == "call") {
        if (a.coherent.empty() || a.out.empty()) { usage(); return 2; }
        call(a);
    } else if (a.command == "strip") {
        if (a.input.empty() || a.output.empty() || a.map.empty()) { usage(); return 2; }
        strip(a);
    } else { usage(); return 2; }
    return 0;
}
