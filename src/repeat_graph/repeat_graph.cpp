//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>
#include <iomanip>

#include "../sequence/overlap.h"
#include "../sequence/vertex_index.h"
#include "../common/config.h"
#include "../common/disjoint_set.h"
#include "repeat_graph.h"


namespace
{
	struct pairhash 
	{
	public:
		template <typename T, typename U>
		std::size_t operator()(const std::pair<T, U> &x) const
		{
			return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
		}
	};

	struct Point2d
	{
		Point2d(FastaRecord::Id curId = FastaRecord::ID_NONE, int32_t curPos = 0, 
			  FastaRecord::Id extId = FastaRecord::ID_NONE, int32_t extPos = 0):
			curId(curId), curPos(curPos), extId(extId), extPos(extPos) {}
		
		FastaRecord::Id curId;
		int32_t curPos;
		FastaRecord::Id extId;
		int32_t extPos;
	};

	struct Point1d
	{
		Point1d(FastaRecord::Id seqId = FastaRecord::ID_NONE, int32_t pos = 0):
			seqId(seqId), pos(pos) {}
		
		FastaRecord::Id seqId;
		int32_t pos;
	};
	
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

bool GraphEdge::isTip() const
{
	return nodeLeft->inEdges.empty() || nodeRight->outEdges.empty();
}

void RepeatGraph::build()
{
	//getting overlaps
	VertexIndex asmIndex(_asmSeqs);
	asmIndex.countKmers(1);
	asmIndex.buildIndex(1, Constants::repeatGraphMaxKmer, 
						Constants::repeatGraphKmerSample);

	OverlapDetector asmOverlapper(_asmSeqs, asmIndex, 
								  Constants::maximumJump, 
								  Parameters::get().minimumOverlap,
								  0);
	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs, false);
	asmOverlaps.findAllOverlaps();

	this->getGluepoints(asmOverlaps);
	this->initializeEdges(asmOverlaps);
}


void RepeatGraph::getGluepoints(const OverlapContainer& asmOverlaps)
{
	Logger::get().debug() << "Computing gluepoints";
	typedef SetNode<Point2d> SetPoint2d;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<SetPoint2d*>> endpoints;
	size_t pointId = 0;

	for (auto& seqOvlps : asmOverlaps.getOverlapIndex())
	{
		for (auto& ovlp : seqOvlps.second)
		{
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curBegin,
										  ovlp.extId, ovlp.extBegin)));
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curEnd,
										  ovlp.extId, ovlp.extEnd)));
		}
	}

	for (auto& seqPoints : endpoints)
	{
		std::sort(seqPoints.second.begin(), seqPoints.second.end(),
				  [](const SetPoint2d* p1, const SetPoint2d* p2)
				  {return p1->data.curPos < p2->data.curPos;});

		for (size_t i = 0; i < seqPoints.second.size() - 1; ++i)
		{
			auto* p1 = seqPoints.second[i];
			auto* p2 = seqPoints.second[i + 1];
			if (abs(p1->data.curPos - p2->data.curPos) < _maxSeparation)
			{
				unionSet(p1, p2);
			}

		}
	}

	std::unordered_map<SetPoint2d*, std::vector<SetPoint2d*>> clusters;
	for (auto seqPoints : endpoints)
	{
		for (auto& endpoint : seqPoints.second)
		{
			clusters[findSet(endpoint)].push_back(endpoint);
		}
	}

	typedef SetNode<Point1d> SetPoint1d;
	std::unordered_map<FastaRecord::Id, std::list<SetPoint1d>> tempGluepoints;
	std::unordered_map<SetPoint1d*, SetPoint1d*> complements;

	for (auto& clustEndpoints : clusters)
	{
		FastaRecord::Id clustSeq = clustEndpoints.second.front()->data.curId;
		if (!clustSeq.strand()) continue;	//only for forward strands

		std::vector<int32_t> positions;
		for (auto& ep : clustEndpoints.second) 
		{
			positions.push_back(ep->data.curPos);
		}
		int32_t clusterXpos = median(positions);

		std::vector<Point1d> clusterPoints;
		clusterPoints.emplace_back(clustSeq, clusterXpos);

		std::list<SetPoint2d> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.emplace_back(ep->data);
		}
		
		//projections on other overlaps
		auto startCmp = [] (const OverlapRange& ovlp, int32_t pos)
						{return ovlp.curBegin < pos;};
		auto& allOvlp = asmOverlaps.getOverlapIndex().at(clustSeq);
		auto leftPos = std::lower_bound(allOvlp.begin(), allOvlp.end(), 
										clusterXpos - 50000, startCmp);
		auto rightPos = std::lower_bound(allOvlp.begin(), allOvlp.end(), 
										 clusterXpos, startCmp);

		for (auto& ovlp = leftPos; ovlp != rightPos; ++ovlp)
		{
			if (ovlp->curBegin <= clusterXpos && clusterXpos <= ovlp->curEnd)
			{
				//TODO: projection with k-mers / alignment
				float lengthRatio = (float)ovlp->extRange() / ovlp->curRange();
				int32_t projectedPos = ovlp->extBegin + 
								float(clusterXpos - ovlp->curBegin) * lengthRatio;
				projectedPos = std::max(ovlp->extBegin, 
										std::min(projectedPos, ovlp->extEnd));

				extCoords.emplace_back(Point2d(clustSeq, clusterXpos,
									   		   ovlp->extId, projectedPos));
			}
		}

		//cluster them
		for (auto& p1 : extCoords)
		{
			for (auto& p2 : extCoords)
			{
				if (p1.data.extId == p2.data.extId &&
					abs(p1.data.extPos - p2.data.extPos) < _maxSeparation)
				{
					unionSet(&p1, &p2);
				}
			}
		}
		std::unordered_map<SetPoint2d*, std::vector<SetPoint2d*>> extClusters;
		for (auto& endpoint : extCoords)
		{
			extClusters[findSet(&endpoint)].push_back(&endpoint);
		}

		//now, get coordinates for each cluster
		for (auto& extClust : extClusters)
		{
			std::vector<int32_t> positions;
			for (auto& ep : extClust.second) 
			{
				positions.push_back(ep->data.extPos);
			}
			int32_t clusterYpos = median(positions);


			FastaRecord::Id extSeq = extClust.second.front()->data.extId;
			clusterPoints.emplace_back(extSeq, clusterYpos);
		}

		//merge with the previous clusters
		std::vector<SetPoint1d*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			int32_t seqLen = _asmSeqs.seqLen(clustPt.seqId);
			Point1d complPt(clustPt.seqId.rc(), seqLen - clustPt.pos - 1);

			bool used = false;
			auto& seqGluepoints = tempGluepoints[clustPt.seqId];
			auto& complGluepoints = tempGluepoints[clustPt.seqId.rc()];
			for (auto& glueNode : seqGluepoints)
			{
				if (abs(glueNode.data.pos - clustPt.pos) < _maxSeparation)
				{
					used = true;
					toMerge.push_back(&glueNode);
				}
			}
			if (!used)
			{
				seqGluepoints.emplace_back(clustPt);
				auto fwdPtr = &seqGluepoints.back();
				complGluepoints.emplace_back(complPt);
				auto revPtr = &complGluepoints.back();

				complements[fwdPtr] = revPtr;
				complements[revPtr] = fwdPtr;
				toMerge.push_back(fwdPtr);
			}
		}
		for (size_t i = 0; i < toMerge.size() - 1; ++i)
		{
			unionSet(toMerge[i], toMerge[i + 1]);
			unionSet(complements[toMerge[i]], 
					 complements[toMerge[i + 1]]);
		}
	}

	std::unordered_map<SetPoint1d*, size_t> setToId;
	for (auto& seqGluepoints : tempGluepoints)
	{
		for (auto& gp : seqGluepoints.second)
		{
			if (!setToId.count(findSet(&gp)))
			{
				setToId[findSet(&gp)] = pointId++;
			}
			_gluePoints[gp.data.seqId].emplace_back(setToId[findSet(&gp)],
													gp.data.seqId, gp.data.pos);
		}
	}

	//for (auto& seqPoints : _gluePoints)
	for (auto& seqRec : _asmSeqs.getIndex())
	{
		auto& seqPoints = _gluePoints[seqRec.first];
		std::sort(seqPoints.begin(), seqPoints.end(),
				  [](const GluePoint& pt1, const GluePoint& pt2)
				  {return pt1.position < pt2.position;});

		//flanking points
		seqPoints.emplace(seqPoints.begin(), pointId++, 
						  seqRec.first, 0);
		seqPoints.emplace_back(pointId++, seqRec.first, 
							   _asmSeqs.seqLen(seqRec.first) - 1);
	}

	for (auto& seqEndpoints : endpoints)
	{
		for (auto ep : seqEndpoints.second) delete ep;
	}
}

void RepeatGraph::initializeEdges(const OverlapContainer& asmOverlaps)
{
	Logger::get().debug() << "Initializing edges";

	typedef std::pair<GraphNode*, GraphNode*> NodePair;
	std::unordered_map<NodePair, std::vector<SequenceSegment>, 
					   pairhash> parallelSegments;
	std::unordered_map<NodePair, NodePair, pairhash> complEdges;

	std::unordered_map<size_t, GraphNode*> nodeIndex;
	auto idToNode = [&nodeIndex, this](size_t nodeId)
	{
		if (!nodeIndex.count(nodeId))
		{
			nodeIndex[nodeId] = this->addNode();
		}
		return nodeIndex[nodeId];
	};

	for (auto& seqEdgesPair : _gluePoints)
	{
		if (!seqEdgesPair.first.strand()) continue;
		FastaRecord::Id complId = seqEdgesPair.first.rc();

		if (seqEdgesPair.second.size() != _gluePoints[complId].size())
		{
			throw std::runtime_error("Graph is not symmetric");
		}

		for (size_t i = 0; i < seqEdgesPair.second.size() - 1; ++i)
		{
			GluePoint gpLeft = seqEdgesPair.second[i];
			GluePoint gpRight = seqEdgesPair.second[i + 1];

			size_t complPos = seqEdgesPair.second.size() - i - 2;
			GluePoint complLeft = _gluePoints[complId][complPos];
			GluePoint complRight = _gluePoints[complId][complPos + 1];

			GraphNode* leftNode = idToNode(gpLeft.pointId);
			GraphNode* rightNode = idToNode(gpRight.pointId);
			NodePair fwdPair = std::make_pair(leftNode, rightNode);

			GraphNode* complLeftNode = idToNode(complLeft.pointId);
			GraphNode* complRightNode = idToNode(complRight.pointId);
			NodePair revPair = std::make_pair(complLeftNode, complRightNode);

			int32_t seqLen = _asmSeqs.seqLen(gpLeft.seqId);
			parallelSegments[fwdPair].emplace_back(gpLeft.seqId, seqLen, 
												   gpLeft.position, 
							  					   gpRight.position);
			parallelSegments[revPair]
				.push_back(parallelSegments[fwdPair].back().complement());

			complEdges[fwdPair] = revPair;
			complEdges[revPair] = fwdPair;
		}
	}

	auto segIntersect = [] (const SequenceSegment& s, const OverlapRange& o)
	{
		return std::min(o.curEnd, s.end) - std::max(o.curBegin, s.start);
	};

	std::unordered_set<NodePair, pairhash> usedPairs;
	for (auto& nodePairSeqs : parallelSegments)
	{
		if (usedPairs.count(nodePairSeqs.first)) continue;
		usedPairs.insert(complEdges[nodePairSeqs.first]);

		//cluster segments based on their overlaps
		std::vector<SetNode<SequenceSegment*>*> segmentsClusters;
		for (auto& seg : nodePairSeqs.second) 
		{
			segmentsClusters.push_back(new SetNode<SequenceSegment*>(&seg));
		}
		for (auto& segOne : segmentsClusters)
		{
			for (auto& segTwo : segmentsClusters)
			{
				//TODO: very inefficient, reimplement
				auto& overlaps = asmOverlaps.getOverlapIndex()
												.at(segOne->data->seqId);
				for (auto& ovlp : overlaps)
				{
					if (ovlp.extId != segTwo->data->seqId) continue;

					int32_t intersectOne = segIntersect(*segOne->data, ovlp);
					int32_t intersectTwo = segIntersect(*segTwo->data, 
														ovlp.reverse());
					float rateOne = (float)intersectOne / 
						(segOne->data->end - segOne->data->start);
					float rateTwo = (float)intersectTwo / 
						(segTwo->data->end - segTwo->data->start);

					if (rateOne > 0.5 && rateTwo > 0.5 &&
						abs(intersectOne - intersectTwo) < _maxSeparation)
					{
						unionSet(segOne, segTwo);
						break;
					}
				}
			}
		}
		std::unordered_map<SetNode<SequenceSegment*>*, 
						   std::vector<SequenceSegment*>> edgeClusters;
		for (auto& setNode : segmentsClusters)
		{
			edgeClusters[findSet(setNode)].push_back(setNode->data);
		}
		//

		//add edge foe each cluster
		std::vector<SequenceSegment> usedSegments;
		for (auto& edgeClust : edgeClusters)
		{
			//in case we have complement edges within the node pair
			auto& anySegment = *edgeClust.second.front();
			if (std::find(usedSegments.begin(), usedSegments.end(), anySegment) 
						  != usedSegments.end()) continue;

			GraphNode* leftNode = nodePairSeqs.first.first;
			GraphNode* rightNode = nodePairSeqs.first.second;
			GraphEdge newEdge(leftNode, rightNode, FastaRecord::Id(_nextEdgeId));
			for (auto& seg : edgeClust.second)
			{
				newEdge.seqSegments.push_back(*seg);
				usedSegments.push_back(seg->complement());
			}

			//check if it's self-complmenet
			bool selfComplement = std::find(usedSegments.begin(), 
						usedSegments.end(), anySegment) != usedSegments.end();
			newEdge.selfComplement = selfComplement;

			this->addEdge(std::move(newEdge));
			if (!selfComplement)
			{
				leftNode = complEdges[nodePairSeqs.first].first;
				rightNode = complEdges[nodePairSeqs.first].second;
				GraphEdge* complEdge = this->addEdge(GraphEdge(leftNode, rightNode, 
												FastaRecord::Id(_nextEdgeId + 1)));
				for (auto& seg : edgeClust.second)
				{
					complEdge->seqSegments.push_back(seg->complement());
				}
			}

			_nextEdgeId += 2;
		}

		for (auto& s : segmentsClusters) delete s;
	}

	this->logEdges();
}


void RepeatGraph::logEdges()
{
	typedef std::pair<SequenceSegment*, GraphEdge*> SegEdgePair;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<SegEdgePair>> sequenceEdges;
	for (auto& edge : this->iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			sequenceEdges[segment.seqId].push_back({&segment, edge});
		}
	}
	for (auto& seqEdgesPair : sequenceEdges)
	{
		std::sort(seqEdgesPair.second.begin(), seqEdgesPair.second.end(),
				  [](const SegEdgePair& s1, const SegEdgePair& s2)
				  	{return s1.first->start < s2.first->start;});
	}
	for (auto& seqEdgesPair : sequenceEdges)
	{
		if (!seqEdgesPair.first.strand()) continue;

		for (size_t i = 0; i < seqEdgesPair.second.size(); ++i)
		{
			SequenceSegment* segment = seqEdgesPair.second[i].first;
			GraphEdge* edge = seqEdgesPair.second[i].second;

			std::string unique = edge->seqSegments.size() == 1 ? "*" : " ";
			Logger::get().debug() << unique << "\t" 
								  << edge->edgeId.signedId() << "\t" 
								  << _asmSeqs.seqName(segment->seqId) << "\t"
								  << segment->start << "\t" 
								  << segment->end << "\t"
								  << segment->end - segment->start;
		}
	}
}


GraphPath RepeatGraph::complementPath(const GraphPath& path)
{
	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(_idToEdge.at((*itEdge)->edgeId.rc()));
	}

	return complEdges;
}

GraphEdge* RepeatGraph::complementEdge(GraphEdge* edge)
{
	return _idToEdge.at(edge->edgeId.rc());
}

GraphNode* RepeatGraph::complementNode(GraphNode* node)
{
	if (!node->outEdges.empty())
	{
		return this->complementEdge(node->outEdges.front())->nodeRight;
	}
	else if(!node->inEdges.empty())
	{
		return this->complementEdge(node->inEdges.front())->nodeLeft;
	}
	return nullptr;
}
