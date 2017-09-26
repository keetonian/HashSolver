#include "genome.hpp"

Genome::Genome(){
  // Anything for default?
}

Genome::Genome(const char* filename){
  read_genome(filename, 0);
}

Genome::Genome(const char* filename, int options){
  read_genome(filename, options);
}

Genome::~Genome(){
  for(uint32_t i = 0; i < genome_vector.size(); i++){
    delete genome_vector.at(i);
  }
  genome_vector.clear();
}

void Genome::read_genome(const char* filename, int options){
  std::ifstream file(filename);
  std::string line;
  if(file.is_open()){
    int i = -1;
    while(getline(file, line)) {
      if(line.size() == 0)
	continue;

      if(line[0] == '>'){
	if(i >= options-1 && options != 0)
	  break;

	if(debug)
	  std::cout << line << std::endl;
	
	genome_vector.push_back(new std::vector<char>());
	genome_names.push_back(line);
	i = genome_vector.size() - 1;
	continue; 
      }

      uint32_t j = 0;
      for(; j < line.size(); j++) {
	char cc = std::toupper(line[j]);
	if(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T' && cc != 'N'){
	  if(debug)
	    std::cout << "Changing character " << cc << " to ";

	  cc = change_base(cc);

	  if(debug)
	    std::cout << cc << std::endl;
	}
	genome_vector.at(i)->push_back(cc);
      }
    }
  }
  else{
    std::cerr << "Unable to open genome" << std::endl;
    exit(0);
  }
  file.close();
  // Frequency counting code.
  /*for(uint32_t i = 0; i < genome_vector.at(0)->size()-100; i++){
    //std::cout << genome_vector.at(0)->at(i);
    if(genome_vector.at(0)->at(i) == 'N')
      continue;
    int a = 0;
    int c = 0;
    int g = 0; 
    int t = 0;
    for(uint32_t j = i; j < i+100; j++){
      switch(genome_vector.at(0)->at(j)){
	case 'A': a++; break;
	case 'C': c++; break;
	case 'G': g++; break;
	case 'T': t++; break;
	default: break;
      }
    }
    std::cout << a << '\t' << c << '\t' << g << '\t' << t << std::endl;
  }*/
}

int Genome::error_check(){
  if(this->genome_vector.size() == 0 || this->genome_names.size() == 0 )
    return 1;
  return 0;
}

char* Genome::getstring(uint32_t index){
  if(error_check())
    return 0;
  uint32_t total = 0;
  uint32_t i = 0;
  for(i = 0; i < genome_vector.size(); i++){
    if(index < total + genome_vector.at(i)->size() && index >= total){
      break;
    }
    total += genome_vector.at(i)->size();
  }

  return getstring(i, index-total);
}

char* Genome::getstring(std::string chromosome, uint32_t index){
  if(error_check())
    return 0;

  uint32_t i = get_chromosome_index(chromosome);
  return getstring(i, index);
}

/*
 *
 * Consider making the return statements into exceptions
 * Also figure out how to give a size to the return char*
 * */
char* Genome::getstring(uint32_t chromosome_index, uint32_t index){
  if(error_check())
    return 0;
  if(chromosome_index > genome_vector.size())
    return 0;
  if(index > genome_vector.at(chromosome_index)->size())
    return 0;

  return &genome_vector.at(chromosome_index)->at(index);
}

void Genome::change(uint32_t index, char c){

}

void Genome::change(std::string chromosome, uint32_t index, char c){
  if(error_check())
    return;

  //uint32_t i = get_chromosome_index(chromosome);

}

int Genome::change_all_N(uint32_t consecutive){
  int changed = 0;
  for(uint32_t i = 0; i < genome_vector.size(); i++){
    std::vector<char> * chromosome = genome_vector.at(i);
    for(uint32_t j = 0; j < chromosome->size(); j++){

      // Find N's
      if(chromosome->at(j) == 'N'){
	uint32_t k = j;

	// Find length of string. This does not deal with long strings of N's well.
	while(k < chromosome->size() && chromosome->at(k) == 'N' && k-j != consecutive)
	  k++;

	// Only change N if it is a spurious string of them
	if(k < chromosome->size() && k-j <= consecutive && ((j > 0 && chromosome->at(j-1) != 'N') || j == 0)){
	  for(uint32_t l = j; l < k; l++){
	    chromosome->at(l) = random_base();
	    changed = 1;
	  }
	}
      }
    }
  }

  return changed;
}

void Genome::write_genome(const char* filename, int options){
  std::ofstream fasta(filename, std::ofstream::out);
  for(uint32_t i = 0; i < genome_vector.size(); i++){
    std::vector<char> * chromosome = genome_vector.at(i);
    if(i != 0)
      fasta << "\n";
    fasta << genome_names.at(i) << std::endl;
    for(uint32_t j = 0; j < chromosome->size(); j++){
      // Add bases to fasta file, with only 60 bases on a line.
      if(j != 0 && j % 60 == 0)
	fasta << "\n";
      fasta << chromosome->at(j); 
    }
  }
  fasta.close();
}

uint32_t Genome::get_chromosome_index(std::string chromosome){
  uint32_t i = 0;
  for(i = 0; i < genome_names.size(); i++){
    if(chromosome == genome_names.at(i))
      break;
  }
  return i;
}

std::vector<std::vector<char> * > * Genome::get_genome(){
  return &genome_vector;
}

std::vector<std::string> * Genome::get_genome_names(){
  return &genome_names;
}

/*
 * Returns a random base
 * */
char Genome::random_base(){
  uint32_t r = rand() % 4;
  switch(r){
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }

  // Should never reach this, but included here for the compiler
  return 0;
}

char Genome::change_base(char cc){
  /*
   * Base codes: from zhanglab.ccmb.med.umich.edu/FASTA/
   *
   * U -> Uridine
   * R -> G A (purine)
   * Y -> T C (pyrimidine)
   * K -> G T (keto)
   * M -> A C (amino)
   * S -> G C (strong)
   * W -> A T (weak)
   * B -> G T C
   * D -> G A T
   * H -> A C T
   * V -> G C A
   * N -> A G C T
   * */
  uint64_t rand_int = rand();
  switch(cc){
    case 'R':
      if(rand_int%2 == 0)
	cc = 'G';
      else
	cc = 'A';
      break;
    case 'Y':
      if(rand_int%2 == 0)
	cc = 'T';
      else
	cc = 'C';
      break;
    case 'K':
      if(rand_int%2 == 0)
	cc = 'G';
      else
	cc = 'T';
      break;
    case 'M':
      if(rand_int%2 == 0)
	cc = 'A';
      else
	cc = 'C';
      break;
    case 'S':
      if(rand_int%2 == 0)
	cc = 'G';
      else
	cc = 'C';
      break;
    case 'W':
      if(rand_int%2 == 0)
	cc = 'A';
      else
	cc = 'T';
      break;
    case 'B':
      if(rand_int%3 == 0)
	cc = 'G';
      else if(rand_int%3 == 1)
	cc = 'T';
      else
	cc = 'C';
      break;
    case 'D':
      if(rand_int%3 == 0)
	cc = 'G';
      else if(rand_int%3 == 1)
	cc = 'A';
      else
	cc = 'T';
      break;
    case 'H':
      if(rand_int%3 == 0)
	cc = 'A';
      else if(rand_int%3 == 1)
	cc = 'C';
      else
	cc = 'T';
      break;
    case 'V':
      if(rand_int%3 == 0)
	cc = 'G';
      else if(rand_int%3 == 1)
	cc = 'C';
      else
	cc = 'A';
      break;
    default:
      std::cerr << "Found unknown character while parsing genome: " << cc << std::endl;
  }
  return cc;
}
