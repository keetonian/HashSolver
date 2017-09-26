#include "commandline.hpp"

char* file;
size_t seed_size;
uint32_t replace_n;
uint32_t l2thresh;

int main(int argc, char** argv){
  int a = parseCommands(argc, argv);

  std::cout << file << std::endl;
  std::cout << seed_size << std::endl;
  std::cout << replace_n << std::endl;
  std::cout << l2thresh << std::endl;

  return 0;

}
