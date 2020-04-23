#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <numeric>

#include <iostream>
#include <seqan/align.h>
#include <typeinfo>
#include <set>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace seqan;
using namespace std;


template <typename TChar, typename TConfig>
void generateText(String<TChar, TConfig> & text, unsigned const length)
{
    unsigned alphabetSize = ValueSize<TChar>::VALUE;
    unsigned minChar = MinValue<TChar>::VALUE;

    resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = std::rand() % alphabetSize - minChar;
    }
}


int main(int argc, char const * argv[])
{
    ArgumentParser parser("Test");

    addOption(parser, ArgParseOption("l", "length", "length of the read", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("e", "errors", "Number of allowed Errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("m", "mutations", "Number of mutations done on the read", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("i", "iterations", "Number of iterations testing search", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("HD", "hammingDistance", ""));

    addOption(parser, ArgParseOption("em", "editMuations", "Allow Insertions and Deletions as mutations"));

    addOption(parser, ArgParseOption("mv", "myversion", "Use my version used in yara"));

//     addOption(parser, ArgParseOption("iv", "itv", "Use my version with in text verification"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "\nAll: \t\t\t\t\t\t"  << "\n";

    return 0;
}



