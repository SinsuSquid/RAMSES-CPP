#include "ramses/Config.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <set>

namespace ramses {

static std::set<std::string> warned_keys;

std::string Config::trim(const std::string& s) const {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::string Config::to_lower(const std::string& s) const {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

bool Config::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[Config] Error: Could not open " << filename << std::endl;
        return false;
    }

    std::string line;
    std::string current_block = "";

    while (std::getline(file, line)) {
        size_t comment_pos = line.find('!');
        if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
        line = trim(line);
        if (line.empty()) continue;

        if (line[0] == '&') {
            current_block = to_lower(line.substr(1));
            continue;
        }
        if (line[0] == '/') {
            current_block = "";
            continue;
        }

        size_t eq_pos = line.find('=');
        if (eq_pos != std::string::npos && !current_block.empty()) {
            std::string key = to_lower(trim(line.substr(0, eq_pos)));
            std::string val = trim(line.substr(eq_pos + 1));
            if (!val.empty() && val.back() == ',') val.pop_back();
            if (!val.empty() && (val.front() == '\'' || val.front() == '"')) val = val.substr(1);
            if (!val.empty() && (val.back() == '\'' || val.back() == '"')) val.pop_back();
            blocks_[current_block][key] = val;
            std::cout << "[Config] Loaded " << current_block << ":" << key << " = " << val << std::endl;
        }
    }
    return true;
}

std::string Config::get(const std::string& block, const std::string& key, const std::string& default_val) const {
    auto b_it = blocks_.find(to_lower(block));
    if (b_it != blocks_.end()) {
        auto k_it = b_it->second.find(to_lower(key));
        if (k_it != b_it->second.end()) return k_it->second;
    }
    
    // Key not found, return default_val
    return default_val;
}

int Config::get_int(const std::string& block, const std::string& key, int default_val) const {
    std::string val = get(block, key, "NOT_FOUND");
    if (val == "NOT_FOUND") return default_val;
    if (val.empty()) return default_val;
    // Handle simple arrays like '3*1' by taking the first element
    size_t star_pos = val.find('*');
    if (star_pos != std::string::npos) val = val.substr(star_pos + 1);
    try { return std::stoi(val); } catch (...) { return default_val; }
}

double Config::get_double(const std::string& block, const std::string& key, double default_val) const {
    std::string val = get(block, key, "NOT_FOUND");
    if (val == "NOT_FOUND") return default_val;
    if (val.empty()) return default_val;
    std::replace(val.begin(), val.end(), 'd', 'e');
    std::replace(val.begin(), val.end(), 'D', 'e');
    try { return std::stod(val); } catch (...) { return default_val; }
}

bool Config::get_bool(const std::string& block, const std::string& key, bool default_val) const {
    std::string val = to_lower(get(block, key, default_val ? "true" : "false"));
    if (val == "true" || val == ".true." || val == "t" || val == "1") return true;
    if (val == "false" || val == ".false." || val == "f" || val == "0") return false;
    return default_val;
}

} // namespace ramses
