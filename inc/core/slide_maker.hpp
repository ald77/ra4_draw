#ifndef H_SLIDE_MAKER
#define H_SLIDE_MAKER

#include <fstream>
#include <string>
#include <vector>

class SlideMaker{
public:
  SlideMaker(std::string fname, std::string aspect_ratio = "169");
  ~SlideMaker() = default;

  void AddSlide(std::vector<std::string> pnames, int ncols = -1, std::string title = "");
  void Close();

private:
  std::ofstream outfile_;
  std::string filename_;
};

#endif
