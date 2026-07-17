#include "ramses/utils/Logger.hpp"
#include <spdlog/sinks/stdout_color_sinks.h>

namespace ramses {

void Logger::init() {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    // Optional: Format pattern: [2026-07-17 12:00:00] [info] message
    console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    
    auto logger = std::make_shared<spdlog::logger>("ramses", console_sink);
    spdlog::set_default_logger(logger);
    
    // Set default level
    spdlog::set_level(spdlog::level::info);
}

int Logger::get_rank() {
#ifdef RAMSES_USE_MPI
    int is_initialized;
    MPI_Initialized(&is_initialized);
    if (is_initialized) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }
#endif
    return 0;
}

} // namespace ramses
