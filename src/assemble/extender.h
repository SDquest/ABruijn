//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <deque>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/contig_generator.h"
#include "chimera.h"

class Extender
{
public:
	Extender(const SequenceContainer& readsContainer, 
			 OverlapContainer& ovlpContainer,
			 int coverage, int genomeSize):
		_readsContainer(readsContainer), 
		_ovlpContainer(ovlpContainer),
		_chimDetector(readsContainer, ovlpContainer),
		_coverage(coverage), _genomeSize(genomeSize),
		_progress(genomeSize)
	{}

	void assembleContigs();
	const std::vector<ContigPath>& getContigPaths() const
		{return _contigPaths;}

private:
	struct ExtensionInfo
	{
		ExtensionInfo(): leftTip(false), rightTip(false),
			numSuspicious(0), meanOverlaps(0), stepsToTurn(0),
			assembledLength(0) {}

		std::vector<FastaRecord::Id> reads;
		bool leftTip;
		bool rightTip;
		int numSuspicious;
		int meanOverlaps;
		int stepsToTurn;
		int assembledLength;
	};

	ExtensionInfo extendContig(FastaRecord::Id startingRead);
	int   countRightExtensions(FastaRecord::Id readId) const;
	bool  extendsRight(const OverlapRange& ovlp) const;
	void  convertToContigs();

	const SequenceContainer& _readsContainer;
	OverlapContainer& _ovlpContainer;
	ChimeraDetector   _chimDetector;
	const int 		  _coverage;
	const int 		  _genomeSize;
	ProgressPercent   _progress;

	std::vector<ExtensionInfo> 	_readLists;
	std::vector<ContigPath> 	_contigPaths;
	cuckoohash_map<FastaRecord::Id, size_t>  	_innerReads;
};
