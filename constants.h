#ifndef IA_CORE_CODER_CONSTANTS_H
#define IA_CORE_CODER_CONSTANTS_H

/*****************************************************************************/
/* constant macros                                                           */
/*****************************************************************************/
#define Q0 1
#define Q1 2
#define Q2 4
#define Q3 8
#define Q4 16
#define Q5 32
#define Q6 64
#define Q7 128
#define Q8 256
#define Q9 512
#define Q10 1024
#define Q11 2048
#define Q12 4096
#define Q13 8192
#define Q14 16384
#define Q15 32768
#define Q16 65536
#define Q17 131072
#define Q18 262144
#define Q19 524288
#define Q20 1048576
#define Q21 2097152
#define Q22 4194304
#define Q23 8388608
#define Q24 16777216
#define Q25 33554432
#define Q26 67108864
#define Q27 134217728
#define Q28 268435456
#define Q29 536870912
#define Q30 1073741824
#define Q31 2147483647
#define Q32 4294967296
#define Q35 34359738368
#define Q38 274877906944
#define Q39 549755813887
#define Q40 Q39

#define MAX_64 (WORD64)0x7fffffffffffffffLL
#define MIN_64 (WORD64)0x8000000000000000LL

#define MAX_32 (WORD32)0x7fffffffL
#define MIN_32 (WORD32)0x80000000L

#define MAX_16 (WORD16)0x7fff
#define MIN_16 (WORD16)0x8000

#define MAX_24 (WORD32)0x007fffffL
#define MIN_24 (WORD32)0xFF800000L

#define NULLPTR ((VOID *)0)

#define IT_NULL ((VOID *)0)

#define ADJ_SCALE 11

/*****************************************************************************/
/* function macros                                                           */
/*****************************************************************************/

#endif /* IA_CORE_CODER_CONSTANTS_H */
