/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "header_printer.h"
#include "PsiOutStream.h"
#include "psi4/psi4-dec.h"

#include <cstdio>
#include <sstream>
#include <algorithm>

namespace psi {

HeaderPrinter::HeaderPrinter(const std::string& title, BannerStyle style, int width)
    : title_(title), style_(style), width_(width) {}

HeaderPrinter& HeaderPrinter::subtitle(const std::string& text) {
    subtitle_ = text;
    return *this;
}

HeaderPrinter& HeaderPrinter::add_line(const std::string& text, bool centered) {
    lines_.emplace_back(Line::Type::TEXT, text, centered);
    return *this;
}

HeaderPrinter& HeaderPrinter::add_parameter(const std::string& key, const std::string& value, int key_width) {
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "    %-*s %11s", key_width, (key + ":").c_str(), value.c_str());
    lines_.emplace_back(Line::Type::PARAMETER, std::string(buffer));
    return *this;
}

HeaderPrinter& HeaderPrinter::add_parameter(const std::string& key, int value, int key_width) {
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "    %-*s %11d", key_width, (key + ":").c_str(), value);
    lines_.emplace_back(Line::Type::PARAMETER, std::string(buffer));
    return *this;
}

HeaderPrinter& HeaderPrinter::add_parameter(const std::string& key, long value, int key_width) {
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "    %-*s %11ld", key_width, (key + ":").c_str(), value);
    lines_.emplace_back(Line::Type::PARAMETER, std::string(buffer));
    return *this;
}

HeaderPrinter& HeaderPrinter::add_parameter(const std::string& key, double value,
                                            const std::string& format, int key_width) {
    char buffer[256];
    char value_str[64];
    snprintf(value_str, sizeof(value_str), format.c_str(), value);
    snprintf(buffer, sizeof(buffer), "    %-*s %s", key_width, (key + ":").c_str(), value_str);
    lines_.emplace_back(Line::Type::PARAMETER, std::string(buffer));
    return *this;
}

HeaderPrinter& HeaderPrinter::add_parameter(const std::string& key, bool value, int key_width) {
    return add_parameter(key, value ? "Yes" : "No", key_width);
}

HeaderPrinter& HeaderPrinter::add_authors(const std::vector<std::string>& authors) {
    authors_ = authors;
    return *this;
}

HeaderPrinter& HeaderPrinter::add_separator() {
    lines_.emplace_back(Line::Type::SEPARATOR);
    return *this;
}

HeaderPrinter& HeaderPrinter::add_blank_line() {
    lines_.emplace_back(Line::Type::BLANK);
    return *this;
}

HeaderPrinter& HeaderPrinter::set_custom_banner(const std::string& top, const std::string& bottom) {
    custom_banner_top_ = top;
    custom_banner_bottom_ = bottom.empty() ? top : bottom;
    style_ = BannerStyle::CUSTOM;
    return *this;
}

void HeaderPrinter::print() const {
    print_to(outfile);
}

void HeaderPrinter::print_if(bool condition) const {
    if (condition) {
        print();
    }
}

void HeaderPrinter::print_to(std::shared_ptr<PsiOutStream> out) const {
    // Print top banner
    print_banner_top(out);

    // Print content
    print_content(out);

    // Print bottom banner
    print_banner_bottom(out);

    // Add a blank line after the header
    out->Printf("\n");
}

void HeaderPrinter::print_banner_top(std::shared_ptr<PsiOutStream> out) const {
    switch (style_) {
        case BannerStyle::NONE:
            break;

        case BannerStyle::ARROW:
            out->Printf("\n");
            out->Printf("  ==> %s <==\n\n", title_.c_str());
            break;

        case BannerStyle::BOX: {
            std::string sep(width_, '-');
            out->Printf("\t %s\n", sep.c_str());
            out->Printf("\t %s\n", center_text(title_, width_).c_str());
            if (!subtitle_.empty()) {
                out->Printf("\t %s\n", center_text(subtitle_, width_).c_str());
            }
            if (!authors_.empty()) {
                out->Printf("\t %s\n", center_text("", width_).c_str());
                for (const auto& author : authors_) {
                    out->Printf("\t %s\n", center_text(author, width_).c_str());
                }
            }
            out->Printf("\t %s\n", sep.c_str());
            out->Printf("\n");
            break;
        }

        case BannerStyle::CUSTOM:
            if (!custom_banner_top_.empty()) {
                out->Printf("%s\n", custom_banner_top_.c_str());
            }
            break;
    }
}

void HeaderPrinter::print_banner_bottom(std::shared_ptr<PsiOutStream> out) const {
    switch (style_) {
        case BannerStyle::NONE:
        case BannerStyle::ARROW:
            break;

        case BannerStyle::BOX: {
            std::string sep(width_, '-');
            out->Printf("\t %s\n", sep.c_str());
            break;
        }

        case BannerStyle::CUSTOM:
            if (!custom_banner_bottom_.empty()) {
                out->Printf("%s\n", custom_banner_bottom_.c_str());
            }
            break;
    }
}

void HeaderPrinter::print_content(std::shared_ptr<PsiOutStream> out) const {
    for (const auto& line : lines_) {
        switch (line.type) {
            case Line::Type::TEXT:
                if (line.centered) {
                    if (style_ == BannerStyle::BOX) {
                        out->Printf("\t %s\n", center_text(line.text, width_).c_str());
                    } else {
                        out->Printf("  %s\n", line.text.c_str());
                    }
                } else {
                    if (style_ == BannerStyle::BOX) {
                        out->Printf("\t %s\n", line.text.c_str());
                    } else {
                        out->Printf("  %s\n", line.text.c_str());
                    }
                }
                break;

            case Line::Type::PARAMETER:
                out->Printf("%s\n", line.text.c_str());
                break;

            case Line::Type::SEPARATOR: {
                if (style_ == BannerStyle::BOX) {
                    std::string sep(width_, '-');
                    out->Printf("\t %s\n", sep.c_str());
                } else {
                    out->Printf("\n");
                }
                break;
            }

            case Line::Type::BLANK:
                out->Printf("\n");
                break;
        }
    }
}

std::string HeaderPrinter::center_text(const std::string& text, int width) const {
    if (text.length() >= static_cast<size_t>(width)) {
        return text;
    }
    int padding = (width - text.length()) / 2;
    std::string result(padding, ' ');
    result += text;
    // Add remaining spaces to reach width
    result.resize(width, ' ');
    return result;
}

}  // namespace psi
