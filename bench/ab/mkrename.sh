#!/bin/sh
# bench/ab/mkrename.sh —— 从某个 flashvec.h **自动生成**符号改名表
#
# 用法: mkrename.sh <header> <prefix> <out>
#
# 为什么自动生成而不是手写：注册表（ab_ops.h）只管"我要测哪些 op"，符号隔离是头文件
# 的事实。头文件里每一个 "fv_xxx(" 形态的名字（= 函数名，含 static inline 的内部件）
# 全部改名为 <prefix>xxx，于是 ab_variant.cpp 里直接写 fv_vec2_add_x4(...) 就自动
# 落到本侧命名空间——**新入口一落地就自动被隔离，不用改任何改名表**。
#
# 只抓"后面跟左括号"的名字：类型（fv_vec2）、常量（FV_EPS）不会误伤。
set -u
hdr=$1
pfx=$2
out=$3

{
    printf '/* 自动生成（mkrename.sh %s %s）；不要手工编辑。 */\n' "$pfx" "$hdr"
    printf '#ifndef AB_RENAME_GENERATED_H\n#define AB_RENAME_GENERATED_H\n'
    grep -ohE '\bfv_[A-Za-z0-9_]+[[:space:]]*\(' "$hdr" |
        sed -E 's/[[:space:]]*\($//' |
        sort -u |
        while IFS= read -r n; do
            [ -n "$n" ] || continue
            printf '#define %s %s%s\n' "$n" "$pfx" "${n#fv_}"
        done
    printf '#endif\n'
} > "$out.tmp" && mv "$out.tmp" "$out"
