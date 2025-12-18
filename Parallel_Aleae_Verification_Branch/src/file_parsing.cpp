#include <iostream>
#include <fstream>
#include <regex>

#include "file_parsing.hpp"

// TOOD: test these functions please

int parse_in_input_file(const std::string &in_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn)
{
    std::ifstream in_file(in_filename);
    if (in_filename.substr(in_filename.length() - 3, 3) != ".in")
    {
        std::cerr << "Error: Input .in file type is invalid for .in file." << std::endl;
        return 1;
    }

    if (!in_file.is_open())
    {
        std::cerr << "Error: Input .in file failed to be opened" << std::endl;
        return 1;
    }

    std::regex init_line_re("([A-Za-z'][A-Za-z0-9.'_]*)[ ]+([0-9]+)[ ]+(LE|LT|GE|GT|N)([ ]+([0-9]*))*");

    std::string temp("");
    std::string chem_name("");
    std::string chem_amount_str("");
    std::string thresh_type_str("");
    std::string thresh_amount_str("");

    unsigned int lineno = 0;
    while (getline(in_file, temp))
    {
        chem_name = std::regex_replace(temp, init_line_re, "$1");
        chem_amount_str = std::regex_replace(temp, init_line_re, "$2");
        thresh_type_str = std::regex_replace(temp, init_line_re, "$3");
        thresh_amount_str = std::regex_replace(temp, init_line_re, "$5");

        if (chem_name == temp || chem_amount_str == temp || thresh_type_str == temp)
        {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Syntax error at line " << lineno + 1 << std::endl;
            return 1;
        }
        else if (thresh_amount_str == "")
        {
            thresh_amount_str = "0";
        }

        unsigned int chem_amount = (unsigned int)std::stoul(chem_amount_str);

        threshold_types thresh_type = THRESH_N;

        if (thresh_type_str == "LT")
        {
            thresh_type = THRESH_LT;
        }
        else if (thresh_type_str == "LE")
        {
            thresh_type = THRESH_LE;
        }
        else if (thresh_type_str == "GE")
        {
            thresh_type = THRESH_GE;
        }
        else if (thresh_type_str == "GT")
        {
            thresh_type = THRESH_GT;
        }

        unsigned int thresh_amount = (unsigned int)std::stoul(thresh_amount_str);

        if (thresh_type == THRESH_N && thresh_amount > 0 || thresh_type < THRESH_N && thresh_amount == 0)
        {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Syntax error at line " << lineno + 1 << std::endl;
            return 1;
        }

        if (chem_str_to_id.count(chem_name) > 0)
        {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Duplicate chemical at line " << lineno + 1 << std::endl;
            return 1;
        }
        else
        {
            chem_str_to_id.insert(std::pair(chem_name, lineno));
        }

        crn.chems.emplace_back(chem_t(chem_name, chem_amount));
        crn.thresholds.emplace_back(thresh_t(thresh_type, thresh_amount));

        lineno++;
    }

    in_file.close();
    return 0;
}

int parse_r_input_file(const std::string &r_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn)
{
    std::ifstream r_file(r_filename);
    if (r_filename.substr(r_filename.length() - 2, 2) != ".r")
    {
        std::cerr << "Error: Input .r file type is invalid for .in file." << std::endl;
        return 1;
    }

    if (!r_file.is_open())
    {
        std::cerr << "Error: Input .r file failed to be opened" << std::endl;
        return 1;
    }

    std::regex integer_re("([1-9][0-9]*)");
    std::regex equation_re("([^:]+)[ ]*:[ ]*([^:]+)[ ]*:[ ]*([0-9]+([.][0-9]+((E|e)[-]{0,1}[0-9]+)*)*)[ \\r]*$");
    std::regex term_re("([A-Za-z'][A-Za-z0-9.'_]*) ([1-9][0-9]*)");

    unsigned int reaction_no = 0;
    std::string reactant_str("");
    std::string product_str("");
    std::string rate_str("");
    std::string temp("");
    while (getline(r_file, temp))
    {
        reactant_str = std::regex_replace(temp, equation_re, "$1");
        product_str = std::regex_replace(temp, equation_re, "$2");
        rate_str = std::regex_replace(temp, equation_re, "$3");

        if (reactant_str == temp || product_str == temp || rate_str == temp)
        {
            std::cerr << "In " << r_filename << std::endl;
            std::cerr << "Syntax error at line " << reaction_no + 1 << std::endl;
            return 1;
        }

        if (rate_str == "0")
        {
            std::cerr << "In " << r_filename << std::endl;
            std::cerr << "Error at line " << reaction_no + 1 << std::endl;
            std::cerr << "Reaction rate must be a positive non-zero value." << std::endl;
            return 1;
        }

        crn.reactions.emplace_back(std::make_unique<reaction_t>());
        crn.reactions[reaction_no]->rate = std::stod(rate_str);

        for (std::sregex_iterator it_reactant = std::sregex_iterator(reactant_str.begin(), reactant_str.end(), term_re);
             it_reactant != std::sregex_iterator(); ++it_reactant)
        {
            std::string cur_term = std::smatch(*it_reactant).str();

            std::string chem_str = std::regex_replace(cur_term, term_re, "$1");
            std::string coeff_st = std::regex_replace(cur_term, term_re, "$2");

            if (!chem_str_to_id.count(chem_str))
            {
                std::cerr << "Missing chemical at line " << reaction_no + 1 << std::endl;
                return 1;
            }

            crn.reactions[reaction_no]->reactants.emplace_back(chem_str_to_id.at(chem_str));
            crn.reactions[reaction_no]->reactant_deltas.emplace_back((unsigned int)std::stoul(coeff_st));
        }

        for (std::sregex_iterator it_product = std::sregex_iterator(product_str.begin(), product_str.end(), term_re);
             it_product != std::sregex_iterator(); ++it_product)
        {
            std::string cur_term = std::smatch(*it_product).str();

            std::string chem_str = std::regex_replace(cur_term, term_re, "$1");
            std::string coeff_st = std::regex_replace(cur_term, term_re, "$2");

            if (!chem_str_to_id.count(chem_str))
            {
                std::cerr << "In " << r_filename << std::endl;
                std::cerr << "Missing chemical at line " << reaction_no + 1 << std::endl;
                return 1;
            }

            crn.reactions[reaction_no]->products.emplace_back(chem_str_to_id.at(chem_str));
            crn.reactions[reaction_no]->product_deltas.emplace_back((unsigned int)std::stoul(coeff_st));
        }

        reaction_no++;
    }

    r_file.close();
    return 0;
}
