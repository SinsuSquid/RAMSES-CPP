#ifndef RAMSES_LOGGER_HPP
#define RAMSES_LOGGER_HPP

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

class Logger {
public:
    static void init();
    static int get_rank();
};

} // namespace ramses

// Standard logging macros (rank 0 only)
#define RAMSES_TRACE(...) do { if (ramses::Logger::get_rank() == 0) { spdlog::trace(__VA_ARGS__); } } while(0)
#define RAMSES_DEBUG(...) do { if (ramses::Logger::get_rank() == 0) { spdlog::debug(__VA_ARGS__); } } while(0)
#define RAMSES_INFO(...)  do { if (ramses::Logger::get_rank() == 0) { spdlog::info(__VA_ARGS__);  } } while(0)
#define RAMSES_WARN(...)  do { if (ramses::Logger::get_rank() == 0) { spdlog::warn(__VA_ARGS__);  } } while(0)
#define RAMSES_ERROR(...) do { if (ramses::Logger::get_rank() == 0) { spdlog::error(__VA_ARGS__); } } while(0)
#define RAMSES_CRITICAL(...) do { if (ramses::Logger::get_rank() == 0) { spdlog::critical(__VA_ARGS__); } } while(0)

// Logging macros for all ranks
#define RAMSES_TRACE_ALL(...) do { spdlog::trace("[Rank {}] {}", ramses::Logger::get_rank(), fmt::format(__VA_ARGS__)); } while(0)
#define RAMSES_DEBUG_ALL(...) do { spdlog::debug("[Rank {}] {}", ramses::Logger::get_rank(), fmt::format(__VA_ARGS__)); } while(0)
#define RAMSES_INFO_ALL(...)  do { spdlog::info("[Rank {}] {}", ramses::Logger::get_rank(), fmt::format(__VA_ARGS__));  } while(0)
#define RAMSES_WARN_ALL(...)  do { spdlog::warn("[Rank {}] {}", ramses::Logger::get_rank(), fmt::format(__VA_ARGS__));  } while(0)
#define RAMSES_ERROR_ALL(...) do { spdlog::error("[Rank {}] {}", ramses::Logger::get_rank(), fmt::format(__VA_ARGS__)); } while(0)

#endif // RAMSES_LOGGER_HPP
