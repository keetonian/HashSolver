#include "hashtable.h"
#include "iostream"

using namespace std;

int main(){
  read_table_from_file("seed14.hashtable");
  cout << "read hashtable" << endl;
  read_locations_from_file("seed14.locations");
  cout << "read files" << endl;

  for(uint32_t i = 0; i < table_size; i++){
    cout << i << '\t' << get_frequency(i) << '\t' << get_offset(i);
    if(get_frequency(i) > 0){
      for(uint32_t j = 0; j < get_frequency(i); j++){
	cout << '\t' << get_location(get_offset(i)+j);
      }
    }
    cout << endl;
  }

  return 0;

}
