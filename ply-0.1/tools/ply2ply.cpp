#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include <tr1/functional>

#include <ply.hpp>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

using namespace std::tr1::placeholders;

class ply_to_ply_converter
{
public:
  typedef int format_type;
  enum format {
    same_format,
    ascii_format,
    binary_format,
    binary_big_endian_format,
    binary_little_endian_format
  };
  ply_to_ply_converter(format_type format) : format_(format) {}
  bool convert(std::istream& istream, std::ostream& ostream);
private:
  void info_callback(const std::string& filename, std::size_t line_number, const std::string& message);
  void warning_callback(const std::string& filename, std::size_t line_number, const std::string& message);
  void error_callback(const std::string& filename, std::size_t line_number, const std::string& message);
  void magic_callback();
  void format_callback(ply::format_type format, const std::string& version);
  void element_begin_callback();
  void element_end_callback();
  std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> > element_definition_callback(const std::string& element_name, std::size_t count);
  template <typename ScalarType> void scalar_property_callback(ScalarType scalar);
  template <typename ScalarType> std::tr1::function<void (ScalarType)> scalar_property_definition_callback(const std::string& element_name, const std::string& property_name);
  template <typename SizeType, typename ScalarType> void list_property_begin_callback(SizeType size);
  template <typename SizeType, typename ScalarType> void list_property_element_callback(ScalarType scalar);
  template <typename SizeType, typename ScalarType> void list_property_end_callback();
  template <typename SizeType, typename ScalarType> std::tr1::tuple<std::tr1::function<void (SizeType)>, std::tr1::function<void (ScalarType)>, std::tr1::function<void ()> > list_property_definition_callback(const std::string& element_name, const std::string& property_name);
  void comment_callback(const std::string& comment);
  void obj_info_callback(const std::string& obj_info);
  bool end_header_callback();
  format_type format_;
  ply::format_type input_format_, output_format_;
  bool bol_;
  std::ostream* ostream_;
};

void ply_to_ply_converter::info_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
  std::cerr << filename << ":" << line_number << ": " << "info: " << message << std::endl;
}

void ply_to_ply_converter::warning_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
  std::cerr << filename << ":" << line_number << ": " << "warning: " << message << std::endl;
}

void ply_to_ply_converter::error_callback(const std::string& filename, std::size_t line_number, const std::string& message)
{
  std::cerr << filename << ":" << line_number << ": " << "error: " << message << std::endl;
}

void ply_to_ply_converter::magic_callback()
{
  (*ostream_) << "ply" << "\n";
}

void ply_to_ply_converter::format_callback(ply::format_type format, const std::string& version)
{
  input_format_ = format;

  switch (format_) {
    case same_format:
      output_format_ = input_format_;
      break;
    case ascii_format:
      output_format_ = ply::ascii_format;
      break;
    case binary_format:
      output_format_ = ply::host_byte_order == ply::little_endian_byte_order ? ply::binary_little_endian_format : ply::binary_big_endian_format;
      break;
    case binary_big_endian_format:
      output_format_ = ply::binary_big_endian_format;
      break;
    case binary_little_endian_format:
      output_format_ = ply::binary_little_endian_format;
      break;
  };

  (*ostream_) << "format ";
  switch (output_format_) {
    case ply::ascii_format:
      (*ostream_) << "ascii";
      break;
    case ply::binary_little_endian_format:
      (*ostream_) << "binary_little_endian";
      break;
    case ply::binary_big_endian_format:
      (*ostream_) << "binary_big_endian";
      break;
  }
  (*ostream_) << " " << version << "\n";
}

void ply_to_ply_converter::element_begin_callback()
{
  if (output_format_ == ply::ascii_format) {
    bol_ = true;
  }
}

void ply_to_ply_converter::element_end_callback()
{
  if (output_format_ == ply::ascii_format) {
    (*ostream_) << "\n";
  }
}

std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> > ply_to_ply_converter::element_definition_callback(const std::string& element_name, std::size_t count)
{
  (*ostream_) << "element " << element_name << " " << count << "\n";
  return std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> >(
    std::tr1::bind(&ply_to_ply_converter::element_begin_callback, this),
    std::tr1::bind(&ply_to_ply_converter::element_end_callback, this)
  );
}

template <typename ScalarType>
void ply_to_ply_converter::scalar_property_callback(ScalarType scalar)
{
  if (output_format_ == ply::ascii_format) {
    using namespace ply::io_operators;
    if (bol_) {
      bol_ = false;
      (*ostream_) << scalar;
    }
    else {
      (*ostream_) << " " << scalar;
    }
  }
  else {
    if (((ply::host_byte_order == ply::little_endian_byte_order) && (output_format_ == ply::binary_big_endian_format))
      || ((ply::host_byte_order == ply::big_endian_byte_order) && (output_format_ == ply::binary_little_endian_format))) {
      ply::swap_byte_order(scalar);
    }
    ostream_->write(reinterpret_cast<char*>(&scalar), sizeof(scalar));
  }
}

template <typename ScalarType>
std::tr1::function<void (ScalarType)> ply_to_ply_converter::scalar_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
  (*ostream_) << "property " << ply::type_traits<ScalarType>::old_name() << " " << property_name << "\n";
  return std::tr1::bind(&ply_to_ply_converter::scalar_property_callback<ScalarType>, this, _1);
}

template <typename SizeType, typename ScalarType>
void ply_to_ply_converter::list_property_begin_callback(SizeType size)
{
  if (output_format_ == ply::ascii_format) {
    using namespace ply::io_operators;
    if (bol_) {
      bol_ = false;
      (*ostream_) << size;
    }
    else {
      (*ostream_) << " " << size;
    }
  }
  else {
    if (((ply::host_byte_order == ply::little_endian_byte_order) && (output_format_ == ply::binary_big_endian_format))
      || ((ply::host_byte_order == ply::big_endian_byte_order) && (output_format_ == ply::binary_little_endian_format))) {
      ply::swap_byte_order(size);
    }
    ostream_->write(reinterpret_cast<char*>(&size), sizeof(size));
  }
}

template <typename SizeType, typename ScalarType>
void ply_to_ply_converter::list_property_element_callback(ScalarType scalar)
{
  if (output_format_ == ply::ascii_format) {
    using namespace ply::io_operators;
    (*ostream_) << " " << scalar;
  }
  else {
    if (((ply::host_byte_order == ply::little_endian_byte_order) && (output_format_ == ply::binary_big_endian_format))
      || ((ply::host_byte_order == ply::big_endian_byte_order) && (output_format_ == ply::binary_little_endian_format))) {
      ply::swap_byte_order(scalar);
    }
    ostream_->write(reinterpret_cast<char*>(&scalar), sizeof(scalar));
  }
}

template <typename SizeType, typename ScalarType>
void ply_to_ply_converter::list_property_end_callback()
{
}

template <typename SizeType, typename ScalarType>
std::tr1::tuple<std::tr1::function<void (SizeType)>, std::tr1::function<void (ScalarType)>, std::tr1::function<void ()> > ply_to_ply_converter::list_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
  (*ostream_) << "property list " << ply::type_traits<SizeType>::old_name() << " " << ply::type_traits<ScalarType>::old_name() << " " << property_name << "\n";
  return std::tr1::tuple<std::tr1::function<void (SizeType)>, std::tr1::function<void (ScalarType)>, std::tr1::function<void ()> >(
    std::tr1::bind(&ply_to_ply_converter::list_property_begin_callback<SizeType, ScalarType>, this, _1),
    std::tr1::bind(&ply_to_ply_converter::list_property_element_callback<SizeType, ScalarType>, this, _1),
    std::tr1::bind(&ply_to_ply_converter::list_property_end_callback<SizeType, ScalarType>, this)
  );
}

void ply_to_ply_converter::comment_callback(const std::string& comment)
{
  (*ostream_) << comment << "\n";
}

void ply_to_ply_converter::obj_info_callback(const std::string& obj_info)
{
  (*ostream_) << obj_info << "\n";
}

bool ply_to_ply_converter::end_header_callback()
{
  (*ostream_) << "end_header" << "\n";
  return true;
}

bool ply_to_ply_converter::convert(std::istream& istream, std::ostream& ostream)
{
  ply::ply_parser::flags_type ply_parser_flags = 0;

  ply::ply_parser ply_parser(ply_parser_flags);

  std::string ifilename;

  ply_parser.info_callback(std::tr1::bind(&ply_to_ply_converter::info_callback, this, std::tr1::ref(ifilename), _1, _2));
  ply_parser.warning_callback(std::tr1::bind(&ply_to_ply_converter::warning_callback, this, std::tr1::ref(ifilename), _1, _2));
  ply_parser.error_callback(std::tr1::bind(&ply_to_ply_converter::error_callback, this, std::tr1::ref(ifilename), _1, _2));

  ply_parser.magic_callback(std::tr1::bind(&ply_to_ply_converter::magic_callback, this));
  ply_parser.format_callback(std::tr1::bind(&ply_to_ply_converter::format_callback, this, _1, _2));
  ply_parser.element_definition_callback(std::tr1::bind(&ply_to_ply_converter::element_definition_callback, this, _1, _2));

  ply::ply_parser::scalar_property_definition_callbacks_type scalar_property_definition_callbacks;

  ply::at<ply::int8>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::int8>, this, _1, _2);
  ply::at<ply::int16>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::int16>, this, _1, _2);
  ply::at<ply::int32>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::int32>, this, _1, _2);
  ply::at<ply::uint8>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::uint8>, this, _1, _2);
  ply::at<ply::uint16>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::uint16>, this, _1, _2);
  ply::at<ply::uint32>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::uint32>, this, _1, _2);
  ply::at<ply::float32>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::float32>, this, _1, _2);
  ply::at<ply::float64>(scalar_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::scalar_property_definition_callback<ply::float64>, this, _1, _2);

  ply_parser.scalar_property_definition_callbacks(scalar_property_definition_callbacks);

  ply::ply_parser::list_property_definition_callbacks_type list_property_definition_callbacks;

  ply::at<ply::uint8, ply::int8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::int8>, this, _1, _2);
  ply::at<ply::uint8, ply::int16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::int16>, this, _1, _2);
  ply::at<ply::uint8, ply::int32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::int32>, this, _1, _2);
  ply::at<ply::uint8, ply::uint8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::uint8>, this, _1, _2);
  ply::at<ply::uint8, ply::uint16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::uint16>, this, _1, _2);
  ply::at<ply::uint8, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::uint32>, this, _1, _2);
  ply::at<ply::uint8, ply::float32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::float32>, this, _1, _2);
  ply::at<ply::uint8, ply::float64>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint8, ply::float64>, this, _1, _2);

  ply::at<ply::uint16, ply::int8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::int8>, this, _1, _2);
  ply::at<ply::uint16, ply::int16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::int16>, this, _1, _2);
  ply::at<ply::uint16, ply::int32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::int32>, this, _1, _2);
  ply::at<ply::uint16, ply::uint8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::uint8>, this, _1, _2);
  ply::at<ply::uint16, ply::uint16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::uint16>, this, _1, _2);
  ply::at<ply::uint16, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::uint32>, this, _1, _2);
  ply::at<ply::uint16, ply::float32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::float32>, this, _1, _2);
  ply::at<ply::uint16, ply::float64>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint16, ply::float64>, this, _1, _2);

  ply::at<ply::uint32, ply::int8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::int8>, this, _1, _2);
  ply::at<ply::uint32, ply::int16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::int16>, this, _1, _2);
  ply::at<ply::uint32, ply::int32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::int32>, this, _1, _2);
  ply::at<ply::uint32, ply::uint8>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::uint8>, this, _1, _2);
  ply::at<ply::uint32, ply::uint16>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::uint16>, this, _1, _2);
  ply::at<ply::uint32, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::uint32>, this, _1, _2);
  ply::at<ply::uint32, ply::float32>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::float32>, this, _1, _2);
  ply::at<ply::uint32, ply::float64>(list_property_definition_callbacks) = std::tr1::bind(&ply_to_ply_converter::list_property_definition_callback<ply::uint32, ply::float64>, this, _1, _2);

  ply_parser.list_property_definition_callbacks(list_property_definition_callbacks);

  ply_parser.comment_callback(std::tr1::bind(&ply_to_ply_converter::comment_callback, this, _1));
  ply_parser.obj_info_callback(std::tr1::bind(&ply_to_ply_converter::obj_info_callback, this, _1));
  ply_parser.end_header_callback(std::tr1::bind(&ply_to_ply_converter::end_header_callback, this));

  ostream_ = &ostream;
  return ply_parser.parse(istream);
}

int main(int argc, char* argv[])
{
  ply_to_ply_converter::format_type ply_to_ply_converter_format = ply_to_ply_converter::same_format;

  int argi;
  for (argi = 1; argi < argc; ++argi) {

    if (argv[argi][0] != '-') {
      break;
    }
    if (argv[argi][1] == 0) {
      ++argi;
      break;
    }
    char short_opt, *long_opt, *opt_arg;
    if (argv[argi][1] != '-') {
      short_opt = argv[argi][1];
      opt_arg = &argv[argi][2];
      long_opt = &argv[argi][2];
      while (*long_opt != '\0') {
        ++long_opt;
      }
    }
    else {
      short_opt = 0;
      long_opt = &argv[argi][2];
      opt_arg = long_opt;
      while ((*opt_arg != '=') && (*opt_arg != '\0')) {
        ++opt_arg;
      }
      if (*opt_arg == '=') {
        *opt_arg++ = '\0';
      }
    }

    if ((short_opt == 'h') || (std::strcmp(long_opt, "help") == 0)) {
      std::cout << "Usage: ply2ply [OPTION] [[INFILE] OUTFILE]\n";
      std::cout << "Parse an PLY file.\n";
      std::cout << "\n";
      std::cout << "  -h, --help           display this help and exit\n";
      std::cout << "  -v, --version        output version information and exit\n";
      std::cout << "  -f, --format=FORMAT  set format\n";
      std::cout << "\n";
      std::cout << "FORMAT may be one of the following: ascii, binary, binary_big_endian,\n";
      std::cout << "binary_little_endian.\n";
      std::cout << "If no format is given, the format of INFILE is kept.\n";
      std::cout << "\n";
      std::cout << "With no INFILE/OUTFILE, or when INFILE/OUTFILE is -, read standard input/output.\n";
      std::cout << "\n";
      std::cout << "Report bugs to <" << PACKAGE_BUGREPORT << ">.\n";
      return EXIT_SUCCESS;
    }

    else if ((short_opt == 'v') || (std::strcmp(long_opt, "version") == 0)) {
      std::cout << "ply2ply (" << PACKAGE_NAME << ") " << PACKAGE_VERSION << "\n";
      std::cout << "Copyright (C) 2007 " << PACKAGE_AUTHOR << "\n";
      std::cout << "\n";
      std::cout << "This program is free software; you can redistribute it and/or modify\n";
      std::cout << "it under the terms of the GNU General Public License as published by\n";
      std::cout << "the Free Software Foundation; either version 2 of the License, or\n";
      std::cout << "(at your option) any later version.\n";
      std::cout << "\n";
      std::cout << "This program is distributed in the hope that it will be useful,\n";
      std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
      std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
      std::cout << "GNU General Public License for more details.\n";
      std::cout << "\n";
      std::cout << "You should have received a copy of the GNU General Public License\n";
      std::cout << "along with this program; if not, write to the Free Software\n";
      std::cout << "Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA\n";
      return EXIT_SUCCESS;
    }

    else if ((short_opt == 'f') || (std::strcmp(long_opt, "format") == 0)) {
      if (strcmp(opt_arg, "ascii") == 0) {
        ply_to_ply_converter_format = ply_to_ply_converter::ascii_format;
      }
      else if (strcmp(opt_arg, "binary") == 0) {
        ply_to_ply_converter_format = ply_to_ply_converter::binary_format;
      }
      else if (strcmp(opt_arg, "binary_little_endian") == 0) {
        ply_to_ply_converter_format = ply_to_ply_converter::binary_little_endian_format;
      }
      else if (strcmp(opt_arg, "binary_big_endian") == 0) {
        ply_to_ply_converter_format = ply_to_ply_converter::binary_big_endian_format;
      }
      else {
        std::cerr << "ply2ply: " << "invalid option `" << argv[argi] << "'" << "\n";
        std::cerr << "Try `" << argv[0] << " --help' for more information.\n";
        return EXIT_FAILURE;
      }
    }

    else {
      std::cerr << "ply2ply: " << "invalid option `" << argv[argi] << "'" << "\n";
      std::cerr << "Try `" << argv[0] << " --help' for more information.\n";
      return EXIT_FAILURE;
    }
  }

  int parc = argc - argi;
  char** parv = argv + argi;
  if (parc > 2) {
    std::cerr << "ply2ply: " << "too many parameters" << "\n";
    std::cerr << "Try `" << argv[0] << " --help' for more information.\n";
    return EXIT_FAILURE;
  }

  std::ifstream ifstream;
  const char* ifilename = "";
  if (parc > 0) {
    ifilename = parv[0];
    if (std::strcmp(ifilename, "-") != 0) {
      ifstream.open(ifilename);
      if (!ifstream.is_open()) {
        std::cerr << "ply2ply: " << ifilename << ": " << "no such file or directory" << "\n";
        return EXIT_FAILURE;
      }
    }
  }

  std::ofstream ofstream;
  const char* ofilename = "";
  if (parc > 1) {
    ofilename = parv[1];
    if (std::strcmp(ofilename, "-") != 0) {
      ofstream.open(ofilename);
      if (!ofstream.is_open()) {
        std::cerr << "ply2ply: " << ofilename << ": " << "could not open file" << "\n";
        return EXIT_FAILURE;
      }
    }
  }

  std::istream& istream = ifstream.is_open() ? ifstream : std::cin;
  std::ostream& ostream = ofstream.is_open() ? ofstream : std::cout;

  class ply_to_ply_converter ply_to_ply_converter(ply_to_ply_converter_format);
  return ply_to_ply_converter.convert(istream, ostream);
}
