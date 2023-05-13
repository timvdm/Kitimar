#pragma once

#include <chrono>

namespace Kitimar {

    class Timer
    {
        public:
            using Clock = std::chrono::high_resolution_clock;
            using TimePoint = std::chrono::time_point<Clock>;
            using Duration = std::chrono::duration<double>;

            void start()
            {
                m_start = Clock::now();
            }

            Duration elapsed() const
            {
                auto end = Clock::now();
                return end - m_start;
            }

            Duration pause()
            {
                auto d = elapsed();
                m_total += d;
                m_start = {};
                return d;
            }

            Duration stop()
            {
                auto d = elapsed();
                m_start = {};
                return d;
            }

            void reset()
            {
                m_start = {};
                m_total = {};
            }

            Duration total() const
            {
                return m_total;
            }

        private:
            TimePoint m_start = {};
            Duration m_total = {};
    };


} // namespace Kitimar
