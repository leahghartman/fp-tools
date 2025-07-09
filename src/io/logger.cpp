#include "logger.h"
#include <iostream>
#include <iomanip>

Logger::Logger(const std::string& filepath) {
  out.open(filepath);
  if (!out.is_open()) {
    throw std::runtime_error("Failed to open log file: " + filepath);
  }
}

Logger::~Logger() {
  if (out.is_open()) {
    out.close();
  }
}

void Logger::info(const std::string& message) {
  out << "  [*] " << message << "\n";
}

void Logger::section(const std::string& title) {
  out << "------------------------------------------------------------------------\n";
  out << title << "\n";
  out << "------------------------------------------------------------------------\n";
}

void Logger::header(const std::string& title) {
  out << "========================================================================\n";
  out << "  " << title << "\n";
  out << "========================================================================\n";
}

void Logger::raw(const std::string& text) {
  out << text;
}

void Logger::newline() {
  out << "\n";
}
