#include <benchmark/benchmark.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename TChar, typename TConfig>
void generateText(seqan::String<TChar, TConfig> & text, unsigned const length)
{
    unsigned alphabetSize = seqan::ValueSize<TChar>::VALUE;
    unsigned minChar = seqan::MinValue<TChar>::VALUE;

    seqan::resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = std::rand() % alphabetSize - minChar;
    }
}



template <typename TAlphabet>
void rank(benchmark::State & state)
{
    unsigned const alphabetSize = seqan::ValueSize<TAlphabet>::VALUE;
    size_t windowSize = 64;

    if (state.range(0) <= static_cast<decltype(state.range(0))>(0))
        throw std::invalid_argument{"Text needs to be at least 1 character long."};

//     uint8_t const log_sigma{static_cast<uint8_t>(std::clamp(std::ceil(std::log2(sigma)), 1.0, 64.0))};
    size_t const text_size{static_cast<size_t>(state.range(0))};

    std::mt19937_64 sigma_engine(12345);
    std::uniform_int_distribution<> sigma_dist(0, alphabetSize - 1);
    auto sigma_gen = [&sigma_dist, &sigma_engine]() { return sigma_dist(sigma_engine); };

    std::mt19937_64 position_engine(4654561);
    std::uniform_int_distribution<> position_dist(0, text_size - 1 - windowSize);
    auto position_gen = [&position_dist, &position_engine]() { return position_dist(position_engine); };

    String<TAlphabet, Packed<> > genome;
    generateText(genome, text_size);

    for (auto _ : state)
    {
        size_t smaller = 0;
        size_t pos = position_gen();
        uint32_t rank;
        TAlphabet val = sigma_gen();

        benchmark::DoNotOptimize(infix(genome, pos, pos + windowSize));
//         benchmark::DoNotOptimize(genome[pos]);

        //         std::cout << "string: " << infix(genome, pos, pos + windowSize) << "\n";
    }



}


BENCHMARK_TEMPLATE(rank, Dna)->RangeMultiplier(100)->Range(100, 100'000'000);
BENCHMARK_TEMPLATE(rank, AminoAcid)->RangeMultiplier(100)->Range(100, 100'000'000);
// BENCHMARK_TEMPLATE(rank, Finite<8>)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, Finite<16>)->RangeMultiplier(100)->Range(100, 1'000'000);
// BENCHMARK_TEMPLATE(rank, Finite<32>)->RangeMultiplier(100)->Range(100, 1'000'000);


/*
BENCHMARK_TEMPLATE(rank, 8)->RangeMultiplier(100)->Range(100, 1'000'000);
BENCHMARK_TEMPLATE(rank, 16)->RangeMultiplier(100)->Range(100, 1'000'000);
BENCHMARK_TEMPLATE(rank, 32)->RangeMultiplier(100)->Range(100, 1'000'000);
BENCHMARK_TEMPLATE(rank, 64)->RangeMultiplier(100)->Range(100, 1'000'000);
BENCHMARK_TEMPLATE(rank, 128)->RangeMultiplier(100)->Range(100, 1'000'000);
BENCHMARK_TEMPLATE(rank, 255)->RangeMultiplier(100)->Range(100, 1'000'000);
*/

BENCHMARK_MAIN();


