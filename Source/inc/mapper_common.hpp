#ifndef MAPPER_COMMON_H_
#define MAPPER_COMMON_H_

enum class SWAFunction { noswa, seqalign, edlib, opal, ssw };
enum class SeedSelection { naive, hobbes, optimal, fasthash };
enum class FilterAlgorithm { none, SHD, MAGNET, QGRAM };

#endif //MAPPER_COMMON_H_
