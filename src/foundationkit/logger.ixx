module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <optional>
#include <string>
#include <vector>
#endif

export module logger;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <vector>;
import <optional>;
import <fstream>;
#endif

export struct Logger
{
  enum class LogLevel
  {
    DEBUG = 0,
    INFO = 1,
    WARNING = 2,
    ERROR = 3
  };

  std::optional<std::string> fileName;
  std::optional<std::ofstream> stream;
  LogLevel level{LogLevel::DEBUG};
  std::vector<std::pair<LogLevel, std::string>> logs;

  Logger();
  Logger(std::string fileName, LogLevel level);

  void flush();

  void log(LogLevel msgLevel, const std::string& msg);

  void debug(const std::string& msg) { log(LogLevel::DEBUG, msg); };
  void info(const std::string& msg) { log(LogLevel::INFO, msg); };
  void warning(const std::string& msg) { log(LogLevel::WARNING, msg); };
  void error(const std::string& msg) { log(LogLevel::ERROR, msg); };

  Logger& operator+=(const Logger& other);
};