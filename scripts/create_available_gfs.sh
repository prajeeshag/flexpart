#!/usr/bin/env bash
set -eu

path=$1

outfile1=$(pwd)/"AVAILABLE"
outfile=$(pwd)/"__AVAILABLE"
> "$outfile"

# remove trailing slash
path=${path%/}

H=${path##*/}            # 00
ymd=$(basename "$(dirname "$path")")  # 20260317
y=${ymd:0:4}
m=${ymd:4:2}
d=${ymd:6:2}
M=00
S=00
cd $path
for f in gfs.t??z.pgrb2.0p25.f*; do
    step=$(echo "$f" | sed -E 's/.*f([0-9]+).*/\1/')
    valid=$(date -u -d "$y-$m-$d $H:$M:$S ${step} hours" +"%Y%m%d %H%M%S")
    printf "%s      %s      ON DISK\n" "$valid" "$f" >> "$outfile"
done

# sort by time
sort -o "$outfile" "$outfile"

cat <<EOF > "$outfile1"
XXXXXX EMPTY LINES XXXXXXXXX
XXXXXX EMPTY LINES XXXXXXXX
YYYYMMDD HHMMSS   name of the file(up to 80 characters)
EOF
cat "$outfile" >> "$outfile1"
rm "$outfile"
