#pragma once
#include <fstream>
#include <string>
#include <sstream>

class Logger {
public:
    Logger(const std::string& filepath);
    ~Logger();

    void info(const std::string& message);
    void section(const std::string& title);
    void header(const std::string& title);
    void raw(const std::string& text);
    void newline();

private:
    std::ofstream out;
};
