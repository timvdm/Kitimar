#pragma once

#include <Kitimar/Util/Util.hpp>
#include <Kitimar/Util/Timer.hpp>

#include <any>
#include <map>
#include <vector>
#include <iostream>


namespace Kitimar {

    // Implementation
    struct Benchmark
    {
        std::string name;
        std::any result = {};
        Util::Timer::Duration elapsed = {};

        // <name> <result> <elapsed>
        std::string toString(int namePadding = 50, int resultPadding = 10) const
        {
            std::stringstream ss;
            ss << Util::pad(name, namePadding);
            ss << Util::pad(Util::toString(result), resultPadding);
            ss << elapsed.count();
            return ss.str();
        }
    };

    // Input
    struct BenchmarkGroup
    {
        std::string name;
        std::vector<Benchmark> benchmarks;

        std::string toString(std::size_t indent = 4, bool indentName = false) const
        {
            std::stringstream ss;
            // name
            if (indentName)
                ss << std::string(indent, ' ');
            ss << name << '\n';
            // benchmarks
            std::string spaces = std::string(indentName ? 2 * indent: indent, ' ');
            for (auto &benchmark : benchmarks)
                ss << spaces << benchmark.toString() << '\n';
            return ss.str();
        }


    };

    // Task
    struct BenchmarkGroups
    {
        std::string name;
        std::vector<BenchmarkGroup> groups;

        std::string toString(std::size_t indent = 4) const
        {
            std::stringstream ss;
            // name
            ss << name << '\n';
            // groups
            for (auto &group : groups)
                ss << group.toString(indent, true) << '\n';
            return ss.str();
        }

        auto totalElapsed() const
        {
            std::map<std::string, Util::Timer::Duration> total;
            for (auto &group: groups)
                for (auto &benchmark : group.benchmarks)
                    total[benchmark.name] += benchmark.elapsed;
            return total;
        }
    };

    template<typename F>
    auto benchmark(F f, const std::string &name = {})
    {
        Util::Timer timer;
        timer.start();
        if constexpr (std::is_same_v<decltype(f()), void>) {
            f();
            return Benchmark{name, 0, timer.elapsed()};
        } else {
            auto result = f();
            return Benchmark{name, result, timer.elapsed()};
        }
    }



} // namespace Kitimar
