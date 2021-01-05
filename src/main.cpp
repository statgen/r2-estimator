
#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <list>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <getopt.h>

float square(float val)
{
  return val * val;
}

bool is_char_space(char c) { return (bool)std::isspace(c); }

std::list<std::pair<std::string, std::vector<std::size_t>>> parse_groups_file(const std::string& sample_groups_path, const std::vector<std::string>& full_sample_list)
{
  std::list<std::pair<std::string, std::vector<std::size_t>>> ret;

  std::ifstream groups_file(sample_groups_path);

  std::unordered_map<std::string, std::set<std::string>> groups;
  std::string line;
  while (std::getline(groups_file, line))
  {
    auto delim_it = std::find_if(line.begin(), line.end(), is_char_space);
    std::string sample_id(line.begin(), delim_it);
    std::string group_id;
    if (delim_it != line.end())
      group_id.assign(delim_it + 1, line.end());

    groups[group_id].emplace(std::move(sample_id));
  }


  for (auto grp = groups.begin(); grp != groups.end(); ++grp)
  {
    ret.emplace_back();
    ret.back().first = grp->first;
    ret.back().second.reserve(full_sample_list.size());

    std::size_t idx = 0;
    for (auto it = full_sample_list.begin(); it != full_sample_list.end(); ++it,++idx)
    {
      if (grp->second.find(*it) != grp->second.end())
      {
        ret.back().second.emplace_back(idx);
      }
    }
    ret.back().second.shrink_to_fit();
  }

  return ret;
}

std::string to_rounded_string(float val)
{
//  std::ostringstream ss;
//  ss << std::fixed << std::setprecision(4) << (std::round(val * 10000.f) / 10000.f);
//  std::string ret = ss.str();

  std::string ret = std::to_string(val);
  ret.erase(ret.find_last_not_of('0')+1);
  ret.erase(ret.find_last_not_of('.')+1);
  return ret;
}

class prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_ = "/dev/stdout";
  savvy::file::format output_format_ = savvy::file::format::sav2;
  float filter_threshold_ = 0.f;
  bool help_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"filter-threshold", required_argument, 0, 't'},
        {"output", required_argument, 0, 'o'},
        {"output-format", required_argument, 0, 'O'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  savvy::file::format output_format() const { return output_format_; }
  float filter_threshold() const { return filter_threshold_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: r2-estimator [opts ...] <in.{sav,bcf,vcf.gz}> \n";
    os << "\n";
    os << " -h, --help              Print usage\n";
    os << " -o, --output            Path to output file (default: /dev/stdout)\n";
    os << " -O, --output-format     Output file format (sav, bcf, or vcf; default: sav)\n";
    os << " -t, --filter-threshold  List of known sites to compute aggregate stats against\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "ho:O:t:", long_options_.data(), &long_index )) != -1)
    {
      std::string optarg_str = optarg ? optarg : "";
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      case 'o':
        output_path_ = optarg_str;
        break;
      case 'O':
        if (optarg_str == "sav")
          output_format_ = savvy::file::format::sav2;
        else if (optarg_str == "bcf")
          output_format_ = savvy::file::format::bcf;
        else if (optarg_str == "vcf")
          output_format_ = savvy::file::format::vcf;
        else
        {
          std::cerr << "Error: invalid --output-format: " << optarg_str << std::endl;
          return false;
        }
        break;
      case 't':
        filter_threshold_ = std::max(0., std::min(1., atof(optarg ? optarg : "")));
        break;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }


  savvy::v2::variant var;
  std::vector<float> hap_dosages;
  savvy::v2::reader input_file(args.input_path());

  if (!input_file)
  {
    std::cerr << "Error: failed to open input file (" << args.input_path() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  std::list<std::pair<std::string, std::vector<std::size_t>>> groups = {}; //parse_groups_file(sample_groups_path, input_file.samples());

  std::vector<std::pair<std::string, std::string>> headers;
  headers.reserve(input_file.headers().size() + 2 + (groups.size() * 2));
  headers.resize(input_file.headers().size());
  bool maf_header_present = false;
  auto next_it = std::copy_if(input_file.headers().begin(), input_file.headers().end(), headers.begin(), [&maf_header_present](const std::pair<std::string, std::string>& h)
  {
    auto details = savvy::parse_header_value(h.second);
    if (details.id == "MAF")
      maf_header_present = true;
    return (details.id != "AF" && details.id  != "R2");
  });

  headers.erase(next_it, headers.end());

  headers.emplace_back("INFO","<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">");
  headers.emplace_back("INFO","<ID=R2,Number=1,Type=Float,Description=\"R-squared Estimate\">");

  for (auto grp = groups.begin(); grp != groups.end(); ++grp)
  {
    headers.emplace_back("INFO","<ID=AF_" + grp->first + ",Number=1,Type=Float,Description=\"Allele Frequency (" + grp->first + ")\">");
    headers.emplace_back("INFO","<ID=R2_" + grp->first + ",Number=1,Type=Float,Description=\"R-squared Estimate (" + grp->first + ")\">");
  }

  savvy::v2::writer output_file(args.output_path(), args.output_format(), headers, input_file.samples());
  if (!output_file)
  {
    std::cerr << "Error: failed to open output file (" << args.output_path() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  while (input_file >> var)
  {
    if (!var.get_format("HDS", hap_dosages))
    {
      std::cerr << "Error: HDS not present\n";
      return EXIT_FAILURE;
    }
    // if you have 'K' haplotypes (not samples), each with a HDS value, then the estimated R2 = Var(K)/[p(1-p)] where p is the allele frequency.

    const std::size_t ploidy = hap_dosages.size() / input_file.samples().size();

    {
      double af = std::accumulate(hap_dosages.begin(), hap_dosages.end(), 0.) / hap_dosages.size();

      double sos = 0.f;

      for (auto it = hap_dosages.begin(); it != hap_dosages.end(); it++)
      {
        sos += square(*it - af);
      }

      double r2 = 0;

      if (af > 0.f && af < 1.f)
        r2 = (sos / hap_dosages.size()) / (af * (1.f - af));

      if (r2 < args.filter_threshold())
        continue;

      var.set_info("AF", (float)af);
      var.set_info("R2", (float)r2);
      if (maf_header_present)
        var.set_info("MAF", static_cast<float>(af > 0.5 ? 1. - af : af));;
    }

    for (auto grp = groups.begin(); grp != groups.end(); ++grp)
    {
      if (grp->second.size())
      {
        double af = 0.f;

        for (auto it = grp->second.begin(); it != grp->second.end(); ++it)
        {
          for (std::size_t h = 0; h < ploidy; ++h)
            af += hap_dosages[*it + h];
        }

        af = af / (grp->second.size() * ploidy);

        double sos = 0.f;

        for (auto it = grp->second.begin(); it != grp->second.end(); it++)
        {
          for (std::size_t h = 0; h < ploidy; ++h)
            sos += square(hap_dosages[*it + h] - af);
        }

        double r2 = 0;

        if (af > 0.f && af < 1.f)
          r2 = (sos / (grp->second.size() * ploidy)) / (af * (1.f - af));

        std::ostringstream af_ss;
        af_ss << af;
        var.set_info("AF_" + grp->first, (float)af);
        var.set_info("R2_" + grp->first, (float)r2);
      }
    }



    output_file << var;


    //std::cerr << std::to_string(r2) << " | " << af_ss.str() << std::endl;

  }

  return output_file.good() && !input_file.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}