//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "multiplicity_inferer.h"
#include "../common/disjoint_set.h"

#include <simplex.h>
#include <variable.h>

void MultiplicityInferer::
	fixEdgesMultiplicity(const std::vector<GraphAlignment>& readAln)
{
	this->estimateByCoverage(readAln);
	//this->balanceGraph();
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

	template<typename T>
	T q75(std::vector<T>& vec)
	{
		std::sort(vec.begin(), vec.end());
		//NOTE: there's a bug in libstdc++ nth_element, 
		//that sometimes leads to a segfault
		//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, 
		//				 vec.end());
		return vec[vec.size() * 3 / 4];
	}
}


void MultiplicityInferer::
	estimateByCoverage(const std::vector<GraphAlignment>& readAln)
{
	const int WINDOW = 100;
	const int SHORT_EDGE = Constants::trustedEdgeLength;

	//alternative coverage
	std::unordered_map<GraphEdge*, std::vector<int>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		int numWindows = edge->length() / WINDOW;
		wndCoverage[edge].assign(numWindows, 0);
	}

	for (auto& path : readAln)
	{
		for (size_t i = 0; i < path.size(); ++i)
		{
			auto& ovlp = path[i].overlap;
			auto& coverage = wndCoverage[path[i].edge];
			for (int pos = ovlp.extBegin / WINDOW + 1; 
			 	 pos < ovlp.extEnd / WINDOW; ++pos)
			{
				if (pos >= 0 && 
					pos < (int)coverage.size())
				{
					++coverage[pos];
				}
			}
		}
	}

	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto& edgeCoverage : wndCoverage)
	{
		if (edgeCoverage.first->length() < SHORT_EDGE) continue;
		for (auto& cov : edgeCoverage.second)
		{
			sumCov += cov;
			++sumLength;
		}
	}
	_meanCoverage = (sumLength != 0) ? sumCov / sumLength : 1;

	Logger::get().debug() << "Mean edge coverage: " << _meanCoverage;

	std::vector<int> edgesCoverage;
	for (auto edge : _graph.iterEdges())
	{
		if (wndCoverage[edge].empty()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		int medianCov = (median(wndCoverage[edge]) + 
						 median(wndCoverage[complEdge])) / 2;

		float minMult = (!edge->isTip()) ? 1 : 0;
		int estMult = std::max(minMult, 
							   roundf((float)medianCov / _meanCoverage));
		if (estMult == 1)
		{
			edgesCoverage.push_back(medianCov);
		}

		std::string match = estMult != edge->multiplicity ? "*" : " ";
		std::string covStr;
		/*for (int cov : altCoverage[edge])
		{
			covStr += std::to_string(cov) + " ";
		}*/
		Logger::get().debug() << match << "\t" << edge->edgeId.signedId() << "\t"
				<< edge->length() << "\t"
				<< edge->multiplicity << "\t" << estMult << "\t"
				<< medianCov << "\t"
				<< (float)medianCov / _meanCoverage;
		//Logger::get().debug() << covStr;

		edge->multiplicity = estMult;
		edge->meanCoverage = medianCov;
	}

	_uniqueCovThreshold = q75(edgesCoverage);
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

/*
void MultiplicityInferer::balanceGraph()
{
	auto trustedEdge = [](const GraphEdge* edge)
	{
		return edge->length() > Constants::trustedEdgeLength && 
			   edge->multiplicity == 1;
	};

	using namespace optimization;

	Logger::get().info() << "Updating edges multiplicity";

	//enumerating edges
	std::unordered_map<GraphEdge*, size_t> edgeToId;
	std::unordered_map<GraphEdge*, SetNode<GraphEdge*>*> edgeClusters;
	std::map<size_t, GraphEdge*> idToEdge;
	size_t numberEdges = 0;
	size_t numTrusted = 0;
	for (auto edge : _graph.iterEdges())
	{
		if (edge->isLooped()) continue;

		if (trustedEdge(edge))
		{
			++numTrusted;
			continue;
		}

		edgeClusters[edge] = new SetNode<GraphEdge*>(edge);

		if (!edgeToId.count(edge))
		{
			GraphEdge* complEdge = _graph.complementPath({edge}).front();
			edgeToId[edge] = numberEdges;
			edgeToId[complEdge] = numberEdges;
			idToEdge[numberEdges] = edge;
			++numberEdges;
		}
	}

	//enumerating nodes
	std::unordered_map<GraphNode*, size_t> nodeToId;
	std::map<size_t, GraphNode*> idToNode;
	size_t numberNodes = 0;
	for (auto node : _graph.iterNodes())
	{
		if (node->inEdges.empty() || node->outEdges.empty()) continue;
		if (node->neighbors().size() < 2) continue;

		for (auto& inEdge : node->inEdges)
		{
			if (trustedEdge(inEdge) || inEdge->isLooped()) continue;
			for (auto& outEdge : node->outEdges)
			{
				if (trustedEdge(outEdge) || outEdge->isLooped()) continue;

				unionSet(edgeClusters[inEdge], 
						 edgeClusters[outEdge]);
			}
		}

		if (!nodeToId.count(node))
		{
			GraphNode* complNode = _graph.complementNode(node);
			nodeToId[complNode] = numberNodes;
			nodeToId[node] = numberNodes;
			idToNode[numberNodes] = node;
			++numberNodes;
		}
	}

	Logger::get().debug() << "Processing graph with " << numberNodes 
		<< " nodes, " << numberEdges << " edges";
	Logger::get().debug() << "Found " << numTrusted / 2 << " trusted edges";

	std::unordered_map<SetNode<GraphEdge*>*, int> clusters;
	for (auto nodePair : edgeClusters) clusters[findSet(nodePair.second)] += 1;
	Logger::get().debug() << "Clusters: " << clusters.size();
	for (auto clSize : clusters)
	{
		Logger::get().debug() << "\t" << clSize.second;
	}

	//formulate linear programming
	Simplex simplex("");
	size_t numVariables = numberEdges + numberNodes * 2;
	for (auto& idEdgePair : idToEdge)
	{
		pilal::Matrix eye(1, numVariables, 0);
		eye(idEdgePair.first) = 1;

		std::string varName = 
				std::to_string(idEdgePair.second->edgeId.signedId());
		simplex.add_variable(new Variable(&simplex, varName.c_str()));
	
		simplex.add_constraint(Constraint(eye, CT_MORE_EQUAL, 
									  (float)idEdgePair.second->multiplicity));
	}

	std::vector<std::vector<int>> incorporatedEquations;
	for (auto& idNodePair : idToNode)
	{
		size_t sourceId = numberEdges + idNodePair.first * 2;
		size_t sinkId = numberEdges + idNodePair.first * 2 + 1;

		//emergency source
		std::string sourceName = std::to_string(idNodePair.first) + "_source";
		simplex.add_variable(new Variable(&simplex, sourceName.c_str()));
		pilal::Matrix sourceMat(1, numVariables, 0);
		sourceMat(sourceId) = 1;
		simplex.add_constraint(Constraint(sourceMat, CT_MORE_EQUAL, 0.0f));

		//emergency sink
		std::string sinkName = std::to_string(idNodePair.first) + "_sink";
		simplex.add_variable(new Variable(&simplex, sinkName.c_str()));
		pilal::Matrix sinkMat(1, numVariables, 0);
		sinkMat(sinkId) = 1;
		simplex.add_constraint(Constraint(sinkMat, CT_MORE_EQUAL, 0.0f));
		
		//adding in/out edges into the problem
		std::vector<int> coefficients(numberEdges, 0);
		int degreeSum = 0;
		for (auto edge : idNodePair.second->inEdges) 
		{
			if (!edge->isLooped())
			{
				if (trustedEdge(edge)) 
				{
					degreeSum -= 1;
				}
				else
				{
					coefficients[edgeToId[edge]] += 1;
				}
			}
		}
		for (auto edge : idNodePair.second->outEdges) 
		{
			if (!edge->isLooped())
			{
				if (trustedEdge(edge)) 
				{
					degreeSum += 1;
				}
				else
				{
					coefficients[edgeToId[edge]] -= 1;
				}
			}
		}
		//

		//build the matrix with all equations and check if it's linearly independend
		pilal::Matrix problemMatrix(incorporatedEquations.size() + 1, 
									numberEdges, 0);
		for (size_t column = 0; column < incorporatedEquations.size(); ++column)
		{
			for (size_t row = 0; row < numberEdges; ++row)
			{
				problemMatrix(column, row) = incorporatedEquations[column][row];
			}
		}
		for (size_t row = 0; row < numberEdges; ++row)
		{
			problemMatrix(incorporatedEquations.size(), row) = coefficients[row];
		}
		if (!problemMatrix.rows_linearly_independent()) continue;
		//

		pilal::Matrix coefMatrix(1, numVariables, 0);
		for (size_t i = 0; i < numberEdges; ++i) 
		{
			coefMatrix(i) = coefficients[i];
		}
		coefMatrix(sourceId) = 1;
		coefMatrix(sinkId) = -1;
		simplex.add_constraint(Constraint(coefMatrix, CT_EQUAL, (float)degreeSum));
		incorporatedEquations.push_back(std::move(coefficients));
	}

    pilal::Matrix costs(1, numVariables, 1.0f);
	for (size_t i = numberEdges; i < numVariables; ++i) costs(i) = 1000.0f;
    simplex.set_objective_function(ObjectiveFunction(OFT_MINIMIZE, costs));   

	simplex.solve();
	if (!simplex.has_solutions() || simplex.must_be_fixed() || 
		simplex.is_unlimited()) throw std::runtime_error("Error while solving LP");

	//simplex.print_solution();
	for (auto& edgeIdPair : edgeToId)
	{
		int inferredMult = simplex.get_solution()(edgeIdPair.second);
		if (edgeIdPair.first->multiplicity != inferredMult)
		{
			Logger::get().debug() << "Mult " 	
				<< edgeIdPair.first->edgeId.signedId() << " " <<
				edgeIdPair.first->multiplicity << " -> " << inferredMult;

			edgeIdPair.first->multiplicity = inferredMult;
		}
	}

	//show warning if the graph remained unbalanced
	int nodesAffected = 0;
	int sumSource = 0;
	int sumSink = 0;
	for (size_t i = 0; i < numberNodes; ++i)
	{
		int nodeSource = simplex.get_solution()(numberEdges + i * 2);
		int nodeSink = simplex.get_solution()(numberEdges + i * 2 + 1);
		sumSource += nodeSource;
		sumSink += nodeSink;
		if (nodeSource + nodeSink > 0)
		{
			++nodesAffected;
		}
	}

	if (nodesAffected)
	{
		Logger::get().warning() << "Could not balance assembly graph in full: "
			<< nodesAffected << " nodes remained, extra source: " << sumSource
			<< " extra sink: " << sumSink;
	}
}*/

