#include <string>
#include <iostream>
using namespace std;

std::string getDNAShape(std::string fastaFilePath, std::string shapeType);

int main(int argc, char **argv) {
  std::string fastaFile(argv[1]);
  std::string shape("MGW");
  if (argc > 2) {
    shape = argv[2];
  }
  cout << getDNAShape(fastaFile, shape) << endl;
}
