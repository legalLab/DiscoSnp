Ttot="$(date +%s)"
echo "Time: ${Ttot}"
echo "bob"
Ttot="$(($(date +%s)-Ttot))"
echo "DiscoSnpRad total time in seconds: ${Ttot}"

