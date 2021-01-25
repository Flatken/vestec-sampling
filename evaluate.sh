cat $1 | grep ":Constructor] Elapsed time in milliseconds :" | cut -d ' ' -f 8 | awk '{ sum += $1 } END { if (NR > 0) print "constructor " sum / NR }'
cat $1 | grep "identify_critical_points] Elapsed time in milliseconds :" | cut -d ' ' -f 8 | awk '{ sum += $1 } END { if (NR > 0) print "Detect " sum / NR }'
cat $1 | grep "reduce" | cut -d ' ' -f 8 | awk '{ sum += $1 } END { if (NR > 0) print "Reduce " sum / NR }'

cat $1 | grep "cleanupDataSet] Elapsed time" | cut -d ' ' -f 8 | awk '{ sum += $1 } END { if (NR > 0) print "clean " sum / NR }'
