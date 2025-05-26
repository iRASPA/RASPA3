/*

 Copyright (c) 2023-2024 Oliver Lau <oliver.lau@gmail.com>

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the “Software”), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or
 sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
 AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.

*/


module;

#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

export module getopt;

export namespace argparser
{

    class argument_required_exception : public std::exception
    {
        std::string exception_message;

    public:
        explicit argument_required_exception(std::string const &option_name)
        {
            exception_message = "Option `" + option_name + "` requires an argument, but none is given.";
        }
        const char *what() const throw()
        {
            return exception_message.c_str();
        }
    };

    class unknown_option_exception : public std::exception
    {
        std::string exception_message;

    public:
        explicit unknown_option_exception(std::string const &option_name)
        {
            exception_message = "Option `" + option_name + "` is unknown.";
        }
        const char *what() const throw()
        {
            return exception_message.c_str();
        }
    };

    class help_requested_exception : public std::exception
    {
        std::string exception_message{};

    public:
        explicit help_requested_exception(){};
        const char *what() const throw()
        {
            return exception_message.c_str();
        }
    };

    namespace util
    {
        template <class T>
        inline void hash_combine(std::size_t &seed, const T &v)
        {
            std::hash<T> hasher;
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        struct options_hash
        {
            std::size_t operator()(std::vector<std::string> const &options) const noexcept
            {
                std::size_t seed = 0;
                for (std::string const &s : options)
                {
                    hash_combine(seed, s);
                }
                return seed;
            }
        };

    }

    class argparser
    {
    public:
        enum arg_type_t
        {
            invalid_argument,
            no_argument,
            required_argument,
            optional_argument,
        };
        using callback_t = std::function<void(std::string const &)>;
        struct arg_options
        {
            std::string help{};
            std::string arg_name{};
            callback_t handler{nullptr};
            arg_type_t arg_type{invalid_argument};

            operator bool() const
            {
                return handler != nullptr;
            }
        };
        struct positional
        {
            std::string help{};
            std::string arg_name{};
            callback_t handler{nullptr};
        };

        argparser() = delete;
        argparser(int argc, char *argv[])
        {
            const std::size_t n = static_cast<std::size_t>(argc);
            args_.reserve(n);
            for (std::size_t i = 1; i < n; ++i)
            {
                args_.push_back(argv[i]);
            }
        }

        argparser &reg(std::vector<std::string> const &options, std::string const &arg_name, arg_type_t arg_type, std::string const &help, callback_t handler)
        {
            options_.emplace(std::make_pair(
                options,
                arg_options{help, arg_name, handler, arg_type}));
            return *this;
        }

        argparser &reg(std::vector<std::string> const &options, arg_type_t arg_type, std::string const &help, callback_t handler)
        {
            return reg(options, std::string{}, arg_type, help, handler);
        }

        [[deprecated]]
        argparser &reg(std::vector<std::string> const &options, arg_type_t arg_type, callback_t handler)
        {
            return reg(options, std::string{}, arg_type, std::string{}, handler);
        }

        argparser &pos(std::string const &arg_name, std::string help, callback_t handler)
        {
            positionals_.emplace_back(positional{help,
                                                 arg_name,
                                                 handler});
            return *this;
        }

        [[deprecated]]
        argparser &pos(std::string const &arg_name, callback_t handler)
        {
            positionals_.emplace_back(positional{std::string{},
                                                 arg_name,
                                                 handler});
            return *this;
        }

        void display_help(void)
        {
            std::cout << info_text_ << "\n"
                      << "Usage:\n\n"
                      << "  " << name_ << ' ';
            if (!options_.empty())
            {
                std::cout << "[OPTIONS]";
                if (!positionals_.empty())
                {
                    std::cout << ' ';
                }
            }
            if (!positionals_.empty())
            {
                bool has_help = false;
                for (auto pos = std::begin(positionals_); pos != std::end(positionals_); ++pos)
                {
                    std::cout << pos->arg_name;
                    if (std::next(pos) != std::end(positionals_))
                    {
                        std::cout << ' ';
                    }
                    has_help |= !pos->help.empty();
                }
                if (has_help)
                {
                    std::cout << "\n\nArguments:\n\n";
                    for (positional const &pos : positionals_)
                    {
                        std::cout << "  " << pos.arg_name << "\n\n"
                                  << "    " << pos.help << "\n";
                    }
                }
            }
            std::cout << "\n\nOptions:\n\n";
            for (auto const &[switches, options] : options_)
            {
                std::cout << "  ";
                for (auto opt = std::begin(switches); opt != std::end(switches); ++opt)
                {
                    switch (options.arg_type)
                    {
                    case no_argument:
                        std::cout << *opt;
                        break;
                    case required_argument:
                        std::cout << *opt << ' ' << options.arg_name;
                        break;
                    case optional_argument:
                        std::cout << *opt << " [" << options.arg_name << "]";
                        break;
                    default:
                        break;
                    }
                    if (std::next(opt) != std::end(switches))
                    {
                        std::cout << ", ";
                    }
                }
                if (!options.help.empty())
                {
                    std::cout << "\n\n    " << options.help;
                }
                std::cout << "\n\n";
            }
            std::cout << std::endl;
            throw help_requested_exception();
        }

        argparser &help(std::vector<std::string> const &options, std::string const &help)
        {
            options_.emplace(std::make_pair(
                options,
                arg_options{help, std::string{}, std::bind(&argparser::display_help, this), no_argument}));
            return *this;
        }

        argparser &info(std::string const &text, std::string const &name)
        {
            info_text_ = text;
            name_ = name;
            return *this;
        }

        arg_options find(std::string const &opt)
        {
            for (auto const &[options, option_options] : options_)
            {
                for (auto const &option : options)
                {
                    if (option == opt)
                    {
                        return option_options;
                    }
                }
            }
            return arg_options{};
        }

        /**
         * @brief Parse command-line arguments.
         *
         * @return true Successfully parsed.
         * @return false Errors occurred.
         */
        bool operator()()
        {
            current_positional_ = positionals_.begin();
            auto arg = args_.begin();
            while (arg != args_.end())
            {
                arg_options opt = find(*arg);
                if (opt)
                {
                    switch (opt.arg_type)
                    {
                    case no_argument:
                        opt.handler(std::string());
                        break;
                    case required_argument:
                    {
                        if (std::next(arg) != std::end(args_))
                        {
                            opt.handler(*(++arg));
                        }
                        else
                        {
                            throw argument_required_exception(*arg);
                        }
                        break;
                    }
                    case optional_argument:
                    {
                        if (std::next(arg) != std::end(args_) && find(*std::next(arg)))
                        {
                            opt.handler(*(++arg));
                        }
                        else
                        {
                            opt.handler(std::string{});
                        }
                        break;
                    }
                    case invalid_argument:
                        [[fallthrough]];
                    default:
                        break;
                    }
                }
                else
                {
                    if (current_positional_ != std::end(positionals_))
                    {
                        (current_positional_->handler)(*arg);
                        ++current_positional_;
                    }
                    else
                    {
                        throw unknown_option_exception(*arg);
                    }
                }
                ++arg;
            }
            return true;
        }

    private:
        std::vector<std::string> args_;
        std::vector<positional> positionals_;
        std::vector<positional>::const_iterator current_positional_;
        std::unordered_map<std::vector<std::string>, arg_options, util::options_hash> options_;
        std::string info_text_;
        std::string name_;
    };

}

