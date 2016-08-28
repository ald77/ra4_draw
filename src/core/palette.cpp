/*! \class Palette

  \brief Loads colors from a text configuration file

  The text files may have multiple palettes defined, with the start of each
  palette marker by a line containing "[palette_name]" and ended by the start of
  the next palette or the end of the file. Within each palette, colors are
  defined on separate lines using their RGB values as "color_name red green
  blue".
*/
#include "core/palette.hpp"

#include <fstream>
#include <sstream>

#include "TColor.h"

#include "core/utilities.hpp"

using namespace std;

/*!\brief Construct from file and palette names

  \param[in] file File from which to read colors

  \param[in] palette %Palette name from which to read colors
*/
Palette::Palette(const string &file,
                 const string &palette):
  file_(file),
  palette_(palette){
  }

/*!\brief Construct from palette name

  Colors are read from the file "txt/colors.txt"

  \param[in] palette %Palette name from which to read colors
*/
Palette::Palette(const string &palette):
  Palette("txt/colors.txt", palette){
}

/*!\brief Get the file from which colors are read

  \return File from which colors are read
*/
const string & Palette::File() const{
  return file_;
}

/*!\brief Set the file from which colors are read

  \param[in] file File from which colors are read

  \return Reference to *this
*/
Palette & Palette::File(const string &file){
  file_ = file;
  return *this;
}

/*!\brief Get the name of the palette from which colors are read

  \return Name of the palette from which colors are read
*/
const string & Palette::PaletteName() const{
  return palette_;
}

/*!\brief Set the name of the palette from which colors are read

  \param[in] palette Name of the palette from which colors are read

  \return Reference to *this
*/
Palette & Palette::PaletteName(const string &palette){
  return PaletteName(file_, palette);
}

/*!\brief Set the file and name of the palette from which colors are read

  \param[in] file File from which colors are read

  \param[in] palette Name of the palette from which colors are read

  \return Reference to *this
*/
Palette & Palette::PaletteName(const string &file,
                               const string &palette){
  file_ = file;
  palette_ = palette;
  return *this;
}

/*!\brief Gets the ROOT color number corresponding to a color in the
  configuration file

  \param[in] color_name Name of color in current palette for which to obtain
  color number

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::operator()(const string &color_name) const{
  ifstream file(file_);
  string line;
  string current_palette_ = "";
  int line_num = 0;
  while(getline(file, line)){
    ++line_num;
    ReplaceAll(line, "=", " ");
    ReplaceAll(line, "\t", " ");
    auto start  = line.find('[');
    auto end = line.find(']');
    if(start==string::npos && end!=string::npos){
      ERROR("Could not find opening brace in line "+to_string(line_num));
    }
    if(start!=string::npos && end==string::npos){
      ERROR("Could not find closing brace in line "+to_string(line_num));
    }
    if(start<end && start != string::npos && end != string::npos){
      current_palette_ = line.substr(start+1, end-start-1);
    }else if(current_palette_ == palette_
             && line.size()
             && line.at(0)!='#'){
      istringstream iss(line);
      string color;
      iss >> color;
      if(color != color_name) continue;
      Int_t r, g, b;
      iss >> r >> g >> b;
      return RGB(r, g, b);
    }
  }
  DBG("No color " << color_name << " in palette " << palette_ << " in file " << file_);
  return 0;
}

/*!\brief Gets the ROOT color number corresponding to a given RGB color

  \param[in] r Red value in range [0, 255]

  \param[in] g Green value in range [0, 255]

  \param[in] b Blue value in range [0, 255]

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::RGB(Int_t r, Int_t g, Int_t b){
  return TColor::GetColor(r, g, b);
}

/*!\brief Gets the ROOT color number corresponding to a given RGB color

  \param[in] r Red value in range [0., 1.]

  \param[in] g Green value in range [0., 1.]

  \param[in] b Blue value in range [0., 1.]

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::RGB(Float_t r, Float_t g, Float_t b){
  return TColor::GetColor(r, g, b);
}

/*!\brief Gets the ROOT color number corresponding to a given HSV color

  \param[in] h Hue in range [0., 360.]

  \param[in] s Saturation in range [0., 1.].

  \param[in] v Value or brightness in range [0., 1.]

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::HSV(Float_t h, Float_t s, Float_t v){
  Float_t r, g, b;
  TColor::HSV2RGB(h, s, v, r, g, b);
  return RGB(r, g, b);
}

/*!\brief Gets the ROOT color number corresponding to a given HLS/HSL color

  \param[in] h Hue in range [0, 255]

  \param[in] l Saturation in range [0, 255].

  \param[in] s Value or brightness in range [0, 255]

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::HLS(Int_t h, Int_t l, Int_t s){
  Int_t r, g, b;
  TColor::HLS2RGB(h, l, s, r, g, b);
  return RGB(r, g, b);
}

/*!\brief Gets the ROOT color number corresponding to a given HLS/HSL color

  \param[in] h Hue in range [0., 360.]

  \param[in] l Saturation in range [0., 1.].

  \param[in] s Value or brightness in range [0., 1.]

  \return Number of color as used by ROOT's TColor system
*/
Int_t Palette::HLS(Float_t h, Float_t l, Float_t s){
  Float_t r, g, b;
  TColor::HLS2RGB(h, l, s, r, g, b);
  return RGB(r, g, b);
}
