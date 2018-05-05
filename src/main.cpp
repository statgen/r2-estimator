
#include <savvy/reader.hpp>
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

int main()
{
  std::string sample_groups_path = "sample_groups.txt";
  std::string input_path = "chr12.dose.vcf.gz";
  std::string output_path = "/dev/stdout";


  savvy::site_info site;
  std::vector<float> hap_dosages;
  savvy::reader input_file(input_path, savvy::fmt::hds);

  std::list<std::pair<std::string, std::vector<std::size_t>>> groups = parse_groups_file(sample_groups_path, input_file.samples());

  std::vector<std::pair<std::string, std::string>> headers(input_file.headers().size() + 2 + (groups.size() * 2));

  headers[0] = {"INFO","<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"};
  headers[1] = {"INFO","<ID=R2,Number=1,Type=Float,Description=\"R-squared Estimate\">"};

  {
    std::size_t i = 0;
    for (auto grp = groups.begin(); grp != groups.end(); ++grp,++i)
    {
      headers[2 + i * 2] = {"INFO","<ID=AF_" + grp->first + ",Number=1,Type=Float,Description=\"Allele Frequency (" + grp->first + ")\">"};
      headers[2 + i * 2 + 1] = {"INFO","<ID=R2_" + grp->first + ",Number=1,Type=Float,Description=\"R-squared Estimate (" + grp->first + ")\">"};
    }
  }

  auto next_it = std::copy_if(input_file.headers().begin(), input_file.headers().end(), headers.begin() + 2 + (groups.size() * 2), [](const std::pair<std::string, std::string>& h)
  {
    auto details = savvy::parse_header_value(h.second);
    return (details.id != "AF" && details.id  != "R2");
  });

  headers.erase(next_it, headers.end());

  savvy::sav::writer output_file(output_path, input_file.samples().begin(), input_file.samples().end(), headers.begin(), headers.end(), savvy::fmt::hds);


  while (input_file.read(site, hap_dosages))
  {
    // if you have 'K' haplotypes (not samples), each with a HDS value, then the estimated R2 = Var(K)/[2p(1-p)] where p is the allele frequency.

    const std::size_t ploidy = hap_dosages.size() / input_file.samples().size();

    std::unordered_map<std::string, std::string> props;
    for (const auto& info_field : input_file.info_fields())
    {
      std::string tmp = site.prop(info_field);
      props[info_field] = site.prop(info_field);
    }

    {
      float af = std::accumulate(hap_dosages.begin(), hap_dosages.end(), 0.f) / hap_dosages.size();

      float sos = 0.f;

      for (auto it = hap_dosages.begin(); it != hap_dosages.end(); it++)
      {
        sos += square(*it - af);
      }

      float r2 = 0;

      if (af > 0.f && af < 1.f)
        r2 = (sos / hap_dosages.size()) / (2.f * af * (1.f - af));

      std::ostringstream af_ss;
      af_ss << af;
      props["AF"] = af_ss.str();
      props["R2"] = std::to_string(r2);
    }

    for (auto grp = groups.begin(); grp != groups.end(); ++grp)
    {
      if (grp->second.size())
      {
        float af = 0.f;

        for (auto it = grp->second.begin(); it != grp->second.end(); ++it)
        {
          for (std::size_t h = 0; h < ploidy; ++h)
            af += hap_dosages[*it + h];
        }

        af = af / (grp->second.size() * ploidy);

        float sos = 0.f;

        for (auto it = grp->second.begin(); it != grp->second.end(); it++)
        {
          for (std::size_t h = 0; h < ploidy; ++h)
            sos += square(hap_dosages[*it + h] - af);
        }

        float r2 = 0;

        if (af > 0.f && af < 1.f)
          r2 = (sos / (grp->second.size() * ploidy)) / (2.f * af * (1.f - af));

        std::ostringstream af_ss;
        af_ss << af;
        props["AF_" + grp->first] = af_ss.str();
        props["R2_" + grp->first] = std::to_string(r2);
      }
    }



    savvy::site_info out_site(std::string(site.chromosome()), site.position(), std::string(site.ref()), std::string(site.alt()), std::move(props));
    output_file.write(out_site, hap_dosages);


    //std::cerr << std::to_string(r2) << " | " << af_ss.str() << std::endl;

  }

  return output_file.good() && !input_file.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}