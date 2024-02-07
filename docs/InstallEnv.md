
## RISC-V THead GCC
```sh
wget https://occ-oss-prod.oss-cn-hangzhou.aliyuncs.com/resource//1698113812618/Xuantie-900-gcc-linux-5.10.4-glibc-x86_64-V2.8.0-20231018.tar.gz

tar -xvf Xuantie-900-gcc-linux-5.10.4-glibc-x86_64-V2.8.0-20231018.tar.gz -C ${HOME}
```
## Qemu
```sh
git clone https://github.com/T-head-Semi/qemu thead_qemu

./configure --prefix=${HOME}/qemu --target-list=riscv64-linux-user --disable-system

make

make install
```