#ifndef RAMSES_CONFIG_HPP
#define RAMSES_CONFIG_HPP

#include <string>
#include <map>
#include <vector>

namespace ramses {

/**
 * @brief Minimal Fortran Namelist-like parser.
 * 
 * Extracts &BLOCK_NAME k1=v1, k2=v2 / structures.
 */
class Config {
public:
    Config() = default;

    /**
     * @brief Parses a RAMSES .nml file.
     */
    bool parse(const std::string& filename);

    /**
     * @brief Retrieves a value from a specific block.
     */
    std::string get(const std::string& block, const std::string& key, const std::string& default_val = "") const;

    int get_int(const std::string& block, const std::string& key, int default_val = 0) const;
    double get_double(const std::string& block, const std::string& key, double default_val = 0.0) const;
    bool get_bool(const std::string& block, const std::string& key, bool default_val = false) const;
    std::vector<double> get_double_array(const std::string& block, const std::string& key) const;
    std::vector<std::string> get_string_array(const std::string& block, const std::string& key) const;

private:
    std::map<std::string, std::map<std::string, std::string>> blocks_;
    
    std::string trim(const std::string& s) const;
    std::string to_lower(const std::string& s) const;
};

} // namespace ramses

#endif // RAMSES_CONFIG_HPP
