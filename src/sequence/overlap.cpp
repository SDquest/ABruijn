//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>
#include <condition_variable>

#include "overlap.h"
#include "../common/config.h"
#include "../common/parallel.h"


//pre-filtering
bool OverlapDetector::goodStart(int32_t curPos, int32_t extPos, 
								int32_t curLen, int32_t extLen) const
{	
	if (_checkOverhang && 
		std::min(curPos, extPos) > _maxOverhang) return false;

	if (extPos > extLen - _minOverlap ||
	    curPos > curLen - _minOverlap) return false;

	return true;
}


OverlapDetector::JumpRes 
OverlapDetector::jumpTest(int32_t curPrev, int32_t curNext,
						  int32_t extPrev, int32_t extNext) const
{
	if (curNext - curPrev > Constants::maximumJump) return J_END;

	if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
		0 < extNext - extPrev && extNext - extPrev < _maxJump)
	{
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maxJump / Constants::closeJumpRate)
		{
			return J_CLOSE;
		}
		if (abs((curNext - curPrev) - (extNext - extPrev)) 
			< _maxJump / Constants::farJumpRate)
		{
			return J_FAR;
		}
	}
	return J_INCONS;
}


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp) const
{

	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	float meanLength = (ovlp.curRange() + ovlp.extRange()) / 2.0f;
	if (lengthDiff > meanLength / Constants::overlapDivergenceRate)
	{
		return false;
	}

	if (_checkOverhang && ovlp.curId != ovlp.extId.rc())
	{
		if (std::min(ovlp.curBegin, ovlp.extBegin) > 
			_maxOverhang) 
		{
			return false;
		}
		if (std::min(ovlp.curLen - ovlp.curEnd, ovlp.extLen - ovlp.extEnd) > 
			_maxOverhang)
		{
			return false;
		}
	}

	return true;
}

namespace
{
	template<typename T>
	T median(std::vector<T>& vec)
	{
		std::sort(vec.begin(), vec.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
		//				 vec.end());
		return vec[vec.size() / 2];
	}
}


void OverlapDetector::
	processKmer(const VertexIndex::ReadPosition& extReadPos,
				std::vector<OverlapDetector::DPRecord>& extPaths,
				int32_t curPos, const FastaRecord& fastaRec) const
 {
	//no trivial matches
	if (extReadPos.readId == fastaRec.id &&
		extReadPos.position == curPos) return;

	size_t extLen = _seqContainer.seqLen(extReadPos.readId);
	size_t curLen = fastaRec.sequence.length();
	int32_t extPos = extReadPos.position;
	if (extLen < (size_t)_minOverlap) return;

	size_t maxCloseId = 0;
	size_t maxFarId = 0;
	int32_t maxCloseLen = 0;
	int32_t maxFarLen = 0;
	bool extendsClose = false;
	bool extendsFar = false;

	//searching for longest possible extension
	for (size_t pathId = 0; pathId < extPaths.size(); ++pathId)
	{
		JumpRes jumpResult = 
			this->jumpTest(extPaths[pathId].ovlp.curEnd, curPos,
						   extPaths[pathId].ovlp.extEnd, extPos);
		int32_t jumpLength = curPos - extPaths[pathId].ovlp.curBegin;

		switch (jumpResult)
		{
			case J_END:
				break;
			case J_INCONS:
				break;
			case J_CLOSE:
				extPaths[pathId].deleted = true;
				if (jumpLength > maxCloseLen)
				{
					extendsClose = true;
					maxCloseId = pathId;	
					maxCloseLen = jumpLength;
				}
				break;
			case J_FAR:
				if (jumpLength > maxFarLen)
				{
					extendsFar = true;
					maxFarId = pathId;
					maxFarLen = jumpLength;
				}
				break;
		}
	}
	//update the best close extension
	if (extendsClose)
	{
		extPaths[maxCloseId].deleted = false;
		extPaths[maxCloseId].ovlp.curEnd = curPos;
		extPaths[maxCloseId].ovlp.extEnd = extPos;
		extPaths[maxCloseId].shifts.push_back(curPos - extPos);
	}
	//update the best far extension, keep the old path as a copy
	if (extendsFar)
	{
		extPaths.push_back(extPaths[maxFarId]);
		extPaths.back().ovlp.curEnd = curPos;
		extPaths.back().ovlp.extEnd = extPos;
		extPaths.back().shifts.push_back(curPos - extPos);
	}
	//if no extensions possible (or there are no active paths), start a new path
	if (!extendsClose && !extendsFar &&
		(this->goodStart(curPos, extPos, curLen, extLen) ||
		fastaRec.id == extReadPos.readId.rc()))	//TODO: temporary bypass overhang
	{
		OverlapRange ovlp(fastaRec.id, extReadPos.readId,
						  curPos, extPos, curLen, extLen);
		extPaths.push_back({ovlp, {}, false});
	}
	//cleaning up
	for (int i = (int)extPaths.size() - 1; i >= 0; --i)
	{
		if (extPaths[i].deleted)
		{
			extPaths[i] = extPaths.back();
			extPaths.pop_back();
		}
	}

 }

std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool uniqueExtensions)
{
	std::unordered_map<FastaRecord::Id, 
					   std::vector<DPRecord>> activePaths;

	//for all kmers in this read
	for (auto curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.isSolid(curKmerPos.kmer)) continue;

		//for all other occurences of this kmer (extension candidates)
		while (_numberRunning > 0) std::this_thread::yield();

		_kmersQueue.clear();
		_pathsQueue.clear();
		_nextJobId = 0;
		_curPos = curKmerPos.position;
		_curFasta = fastaRec;
		FastaRecord::Id prevReadId = FastaRecord::ID_NONE;
		for (auto extKmerPos : _vertexIndex.byKmer(curKmerPos.kmer))
		{
			if (extKmerPos.readId != prevReadId)
			{
				activePaths[extKmerPos.readId];
				_kmersQueue.push_back(extKmerPos);
				_pathsQueue.push_back(&activePaths[extKmerPos.readId]);
				prevReadId = extKmerPos.readId;
			}
		}

		_workersPaused = false;
		while (!_workersPaused) std::this_thread::yield();
	}

	//post-processing
	std::vector<OverlapRange> detectedOverlaps;
	for (auto& ap : activePaths)
	{
		size_t extLen = _seqContainer.seqLen(ap.first);
		OverlapRange maxOverlap;
		bool passedTest = false;
		for (auto& dpRec : ap.second)
		{
			if (this->overlapTest(dpRec.ovlp))
			{
				dpRec.ovlp.leftShift = median(dpRec.shifts);
				dpRec.ovlp.rightShift = extLen - fastaRec.sequence.length() + 
										dpRec.ovlp.leftShift;

				if (!uniqueExtensions)
				{
					detectedOverlaps.push_back(dpRec.ovlp);
				}
				else
				{
					passedTest = true;
					if (dpRec.ovlp.curRange() > maxOverlap.curRange())
					{
						maxOverlap = dpRec.ovlp;
					}
				}
			}
		}
		if (uniqueExtensions && passedTest)
		{
			detectedOverlaps.push_back(maxOverlap);
		}
	}

	return detectedOverlaps;
}

void OverlapDetector::threadWorker()
{
	bool isRunning = false;
	while (true)
	{
		while (_workersPaused) 
		{
			if (isRunning)
			{
				isRunning = false;
				_numberRunning -= 1;
			}
			std::this_thread::yield();
		}
		if (!isRunning)
		{
			isRunning = true;
			_numberRunning += 1;
		}

		if (_threadsExit) return;
		size_t expected = _nextJobId;
		if (expected == _kmersQueue.size()) 
		{
			_workersPaused = true;
			continue;
		}
		if (_nextJobId.compare_exchange_weak(expected, expected + 1))
		{
			processKmer(_kmersQueue[expected], *_pathsQueue[expected], 
						_curPos, _curFasta);
		}
	}
}

const std::vector<OverlapRange>&
OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	if (!_cached.count(readId))
	{
		const FastaRecord& record = _queryContainer.getIndex().at(readId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(record, _onlyMax);
		this->storeOverlaps(overlaps, readId);	
	}
	return _overlapIndex[readId];
}

void OverlapContainer::storeOverlaps(const std::vector<OverlapRange>& overlaps,
									 FastaRecord::Id seqId)
{
	_cached.insert(seqId);
	_cached.insert(seqId.rc());

	auto& fwdOverlaps = _overlapIndex[seqId];
	auto& revOverlaps = _overlapIndex[seqId.rc()];

	std::unordered_set<FastaRecord::Id> extisting;
	if (_onlyMax)
	{
		for (auto& ovlp : fwdOverlaps) extisting.insert(ovlp.extId);
	}

	for (auto& ovlp : overlaps)
	{
		if (_onlyMax && extisting.count(ovlp.extId)) continue;

		auto revOvlp = ovlp.reverse();
		fwdOverlaps.push_back(ovlp);
		revOverlaps.push_back(ovlp.complement());
		_overlapIndex[revOvlp.curId].push_back(revOvlp);
		_overlapIndex[revOvlp.curId.rc()].push_back(revOvlp.complement());
	}
}

void OverlapContainer::findAllOverlaps()
{
	Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (auto& hashPair : _queryContainer.getIndex())
	{
		allQueries.push_back(hashPair.first);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		auto& fastaRec = _queryContainer.getIndex().at(seqId);
		auto overlaps = _ovlpDetect.getSeqOverlaps(fastaRec, false);

		indexMutex.lock();
		this->storeOverlaps(overlaps, seqId);
		indexMutex.unlock();
	};

	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);

	this->filterOverlaps();
}


void OverlapContainer::filterOverlaps()
{
	Logger::get().debug() << "Filtering overlaps";

	//filter identical overlaps
	int filteredIdent = 0;
	for (auto seqId : _queryContainer.getIndex())
	{
		auto& overlaps = _overlapIndex[seqId.first];
		std::vector<OverlapRange> filtered;
		for (auto& ovlp : overlaps)
		{
			bool found = false;
			for (auto& otherOvlp : filtered)
			{
				if (ovlp.equals(otherOvlp))
				{
					found = true;
					break;
				}
			}
			if (!found) filtered.push_back(ovlp);
		}
		filteredIdent += overlaps.size() - filtered.size();
		overlaps = std::move(filtered);
	}

	//filter contained overlaps
	int filteredContained = 0;
	std::unordered_set<OverlapRange*> contained;
	for (auto seqId : _queryContainer.getIndex())
	{
		auto& overlaps = _overlapIndex[seqId.first];
		for (auto& ovlp : overlaps)
		{
			for (auto& otherOvlp : overlaps)
			{
				if (ovlp.containedBy(otherOvlp))
				{
					contained.insert(&ovlp);
					break;
				}
			}
		}
		std::vector<OverlapRange> filtered;
		for (auto& ovlp : _overlapIndex[seqId.first])
		{
			if (!contained.count(&ovlp)) filtered.push_back(ovlp);
		}
		filteredContained += overlaps.size() - filtered.size();
		overlaps = std::move(filtered);
	}

	Logger::get().debug() << "Filtered " << filteredIdent << " identical and "
		<< filteredContained << " contained overlaps";
}


void OverlapContainer::saveOverlaps(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	for (auto& hashPair : _overlapIndex)
	{
		for (auto& ovlp : hashPair.second)
		{
			fout << ovlp.serialize() << std::endl;
		}
	}
}

void OverlapContainer::loadOverlaps(const std::string& filename)
{
	std::ifstream fin(filename);
	if (!fin.is_open())
	{
		throw std::runtime_error("Can't open overlaps file");
	}

	std::string buffer;
	while(!fin.eof())
	{
		std::getline(fin, buffer);
		if (buffer.empty()) break;
		OverlapRange ovlp;
		ovlp.unserialize(buffer);
		_overlapIndex[ovlp.curId].push_back(ovlp);
	}
}
