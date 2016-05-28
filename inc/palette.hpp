#ifndef H_PALETTE
#define H_PALETTE

#include <string>

#include "RtypesCore.h"

class Palette{
public:
  Palette(const std::string &file,
          const std::string &palette);
  explicit Palette(const std::string &palette = "default");
  Palette(const Palette &) = default;
  Palette & operator=(const Palette &) = default;
  Palette(Palette &&) = default;
  Palette & operator=(Palette &&) = default;
  ~Palette() = default;

  const std::string & File() const;
  Palette & File(const std::string &file);

  const std::string & PaletteName() const;
  Palette & PaletteName(const std::string &palette);
  Palette & PaletteName(const std::string &file,
                        const std::string &palette);

  Int_t operator()(const std::string &color_name) const;

  static Int_t RGB(Int_t r, Int_t g, Int_t b);
  static Int_t RGB(Float_t r, Float_t g, Float_t b);

  static Int_t HSV(Float_t h, Float_t s, Float_t v);

  static Int_t HLS(Int_t h, Int_t l, Int_t s);
  static Int_t HLS(Float_t h, Float_t l, Float_t s);

private:
  std::string file_;//!<File from which to read color definitions
  std::string palette_;//!<Palette name from which to read color definitions
};

#endif
