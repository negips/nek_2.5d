#!/bin/sh

genbox << EOF
taylor.box
EOF
#
reatore2 << EOF
box
taylor
EOF
#
rm -f box.rea
rm -f taylor.rea
genmap << EOF
taylor
0.00001
EOF

