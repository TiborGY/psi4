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

#ifndef _psi_src_lib_libpsi4util_header_printer_h_
#define _psi_src_lib_libpsi4util_header_printer_h_

#include <string>
#include <vector>
#include <memory>

namespace psi {

class PsiOutStream;

/**
 * @brief Utility class for printing formatted headers in Psi4 output
 *
 * This class provides a flexible interface for printing headers with various
 * styles and content types, consolidating the many duplicate print_header()
 * implementations across the codebase.
 *
 * Usage example:
 * @code
 * HeaderPrinter header("My Algorithm", HeaderPrinter::BannerStyle::ARROW);
 * header.add_parameter("Threads", 4)
 *       .add_parameter("Memory [MiB]", 1024)
 *       .add_parameter("Cutoff", 1.0e-12, "%11.0E")
 *       .print_if(print_level_ > 0);
 * @endcode
 */
class HeaderPrinter {
public:
    /**
     * @brief Banner style for the header
     */
    enum class BannerStyle {
        /// No banner, just print content
        NONE,
        /// Standard "==> Title <==" format
        ARROW,
        /// Full box with dashes (like DF-MP2)
        BOX,
        /// Custom banner (user-provided)
        CUSTOM
    };

    /**
     * @brief Construct a new Header Printer object
     *
     * @param title The main title for the header
     * @param style The banner style to use (default: ARROW)
     * @param width The width of the header (default: 60)
     */
    explicit HeaderPrinter(const std::string& title,
                          BannerStyle style = BannerStyle::ARROW,
                          int width = 60);

    /**
     * @brief Set a subtitle (appears below title in BOX style)
     *
     * @param text The subtitle text
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& subtitle(const std::string& text);

    /**
     * @brief Add a custom line of text
     *
     * @param text The text to add
     * @param centered Whether to center the text (default: false)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_line(const std::string& text, bool centered = false);

    /**
     * @brief Add a string parameter (key-value pair)
     *
     * @param key The parameter name
     * @param value The parameter value
     * @param key_width Width for the key column (default: 20)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_parameter(const std::string& key, const std::string& value, int key_width = 20);

    /**
     * @brief Add an integer parameter
     *
     * @param key The parameter name
     * @param value The parameter value
     * @param key_width Width for the key column (default: 20)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_parameter(const std::string& key, int value, int key_width = 20);

    /**
     * @brief Add a long integer parameter
     *
     * @param key The parameter name
     * @param value The parameter value
     * @param key_width Width for the key column (default: 20)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_parameter(const std::string& key, long value, int key_width = 20);

    /**
     * @brief Add a double parameter with custom formatting
     *
     * @param key The parameter name
     * @param value The parameter value
     * @param format Printf-style format string (default: "%11.3E")
     * @param key_width Width for the key column (default: 20)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_parameter(const std::string& key, double value,
                                 const std::string& format = "%11.3E", int key_width = 20);

    /**
     * @brief Add a boolean parameter (prints as "Yes"/"No")
     *
     * @param key The parameter name
     * @param value The parameter value
     * @param key_width Width for the key column (default: 20)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_parameter(const std::string& key, bool value, int key_width = 20);

    /**
     * @brief Add author credits (centered text in BOX style)
     *
     * @param authors The author names
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_authors(const std::vector<std::string>& authors);

    /**
     * @brief Add a separator line
     *
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_separator();

    /**
     * @brief Add a blank line
     *
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& add_blank_line();

    /**
     * @brief Set custom top/bottom banners (for CUSTOM style)
     *
     * @param top The top banner line
     * @param bottom The bottom banner line (if empty, uses top)
     * @return HeaderPrinter& Reference to this object for chaining
     */
    HeaderPrinter& set_custom_banner(const std::string& top, const std::string& bottom = "");

    /**
     * @brief Print the header to outfile
     */
    void print() const;

    /**
     * @brief Print the header to outfile if condition is true
     *
     * @param condition The condition to check
     */
    void print_if(bool condition) const;

    /**
     * @brief Print the header to a specific output stream
     *
     * @param out The output stream to print to
     */
    void print_to(std::shared_ptr<PsiOutStream> out) const;

private:
    struct Line {
        enum class Type { TEXT, PARAMETER, SEPARATOR, BLANK };
        Type type;
        std::string text;
        bool centered;

        Line(Type t, const std::string& txt = "", bool center = false)
            : type(t), text(txt), centered(center) {}
    };

    std::string title_;
    std::string subtitle_;
    BannerStyle style_;
    int width_;
    std::vector<Line> lines_;
    std::vector<std::string> authors_;
    std::string custom_banner_top_;
    std::string custom_banner_bottom_;

    void print_banner_top(std::shared_ptr<PsiOutStream> out) const;
    void print_banner_bottom(std::shared_ptr<PsiOutStream> out) const;
    void print_content(std::shared_ptr<PsiOutStream> out) const;
    std::string center_text(const std::string& text, int width) const;
};

}  // namespace psi

#endif  // _psi_src_lib_libpsi4util_header_printer_h_
