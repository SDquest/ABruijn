//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"


void VertexIndex::countKmers(size_t hardThreshold)
{
	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0 || hardThreshold > 100) 
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}

	//TODO: precount size should correlate with kmer size
	const size_t PRE_COUNT_SIZE = 1024 * 1024 * 1024;
	std::vector<unsigned char> preCounters(PRE_COUNT_SIZE, 0);

	//filling up bloom filter
	Logger::get().info() << "Counting kmers (1/2):";
	ProgressPercent bloomProg(_seqContainer.getIndex().size());
	for (auto& seqPair : _seqContainer.getIndex())
	{
		bloomProg.advance();
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(seqPair.first)))
		{
			if (preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE] != 
				std::numeric_limits<unsigned char>::max())
				++preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
		}
	}

	//counting only kmers that have passed the filter
	Logger::get().info() << "Counting kmers (2/2):";

	std::function<void(const FastaRecord::Id&)> countUpdate = 
	[&preCounters, hardThreshold, this] (const FastaRecord::Id& readId)
	{
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			size_t count = preCounters[kmerPos.kmer.hash() % PRE_COUNT_SIZE];
			if (count >= hardThreshold)
			{
				_kmerCounts.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _seqContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}
	processInParallel(allReads, countUpdate, 
					  Parameters::get().numThreads, true);
	
	for (auto kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
	}
}


void VertexIndex::buildIndex(int minCoverage, int maxCoverage, int filterRatio)
{
	Logger::get().info() << "Building kmer index";

	size_t kmerEntries = 0;
	for (auto kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second && 
			kmer.second <= (size_t)maxCoverage)
		{
			kmerEntries += kmer.second;
		}
	}
	Logger::get().debug() << "Kmer index size: " << kmerEntries;

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[minCoverage, maxCoverage, filterRatio, this] 
	(const FastaRecord::Id& readId)
	{
		int counter = 0;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if ((counter++) % filterRatio) continue;

			size_t count = 0;
			_kmerCounts.find(kmerPos.kmer, count);

			//if ((size_t)minCoverage <= count && count <= 10UL * (size_t)maxCoverage)
			if ((size_t)minCoverage <= count && count <= (size_t)maxCoverage)
			{
				_kmerIndex.insert(kmerPos.kmer, nullptr);
				_kmerIndex.update_fn(kmerPos.kmer, 
					[readId, &kmerPos, count](ReadVector*& vec)
					{
						if (vec == nullptr)
						{
							vec = new ReadVector;
							vec->reserve(count);
						}
						vec->emplace_back(readId, kmerPos.position);
					});
			}
			//if (count > (size_t)maxCoverage)
			//{
			//	_repetitiveKmers.insert(kmerPos.kmer);
			//}
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (auto& hashPair : _seqContainer.getIndex())
	{
		allReads.push_back(hashPair.first);
	}
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, true);

	//sorting by read and position
	for (auto kmerHash : _kmerIndex.lock_table())
	{
		std::sort(kmerHash.second->begin(), kmerHash.second->end(),
				  [](const ReadPosition& rp1, const ReadPosition& rp2)
				  	{return (rp1.readId != rp2.readId) ? rp1.readId < rp2.readId : 
					 		rp1.position < rp2.position;});
	}

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}


void VertexIndex::clear()
{
	for (auto kmerHash : _kmerIndex.lock_table())
	{
		delete kmerHash.second;
	}
	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	_repetitiveKmers.clear();
	_repetitiveKmers.reserve(0);
}
