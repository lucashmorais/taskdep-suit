
#!/bin/bash

cd src

make distclean

export LIBS="${LIBS} -lm"

./configure --prefix=/usr --disable-shared

make

make install

make distclean

cd ..
