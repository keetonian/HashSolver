#include "hashtable.h"
#include "iostream"

using namespace std;

int main(){
  read_table_from_file("seed16-1.hashtable");
  cout << "read hashtable" << endl;
  read_locations_from_file("seed16-1.locations");
  cout << "read files" << endl;

  uint64_t total = 0;
  for(uint32_t i = 0; i < table_size; i++){
    //cout << i << '\t' << get_frequency(i) << '\t' << get_offset(i);
    if(get_frequency(i) > 0){
      for(uint32_t j = 0; j < get_frequency(i); j++){
	//cout << '\t' << get_location(get_offset(i)+j);
	total++;
      }
    }
    //cout << endl;
  }

  cout << total << endl;

  return 0;

}
