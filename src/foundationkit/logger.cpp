module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <ios>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#endif

module logger;

#ifndef USE_LEGACY_HEADERS
import <ios>;
import <string>;
import <iostream>;
import <fstream>;
import <sstream>;
import <vector>;
import <optional>;
#endif

Logger::Logger() {}

Logger::Logger(std::string fileName, LogLevel level)
    : fileName(std::move(fileName)),
      stream(std::make_optional<std::ofstream>(this->fileName.value(), std::ios::out)),
      level(level)
{
  if (stream && !(*stream))
  {
    throw std::ios_base::failure("Failed to open log file: " + fileName);
  }
}

void Logger::flush()
{
  for (auto& logEntry : logs)
  {
    if (logEntry.first <= level)
    {
      (*stream) << logEntry.second;
    }
  }
  (*stream) << std::flush;
}

void Logger::log(LogLevel msgLevel, const std::string& msg) { logs.push_back(std::make_pair(msgLevel, msg)); }

Logger& Logger::operator+=(const Logger& other)
{
  for (auto& logEntry : other.logs)
  {
    log(logEntry.first, logEntry.second);
  }
  return *this;
}