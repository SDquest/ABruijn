//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <algorithm>
#include <unordered_set>

#include "homo_polisher.h"
#include "matrix.h"


namespace
{
	//Computes global pairwise alignment with custom substitution matrix
	//and returns both alignment and score
	float pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						    const SubstitutionMatrix& subsMat,
						    std::string& outOne, std::string& outTwo)
	{
		Matrix<float> scoreMat(seqOne.length() + 1, seqTwo.length() + 1);
		Matrix<char> backtrackMat(seqOne.length() + 1, seqTwo.length() + 1);

		scoreMat.at(0, 0) = 0.0f;
		backtrackMat.at(0, 0) = 0;
		for (size_t i = 0; i < seqOne.length(); ++i) 
		{
			scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + 
									subsMat.getScore(seqOne[i], '-');
			backtrackMat.at(i + 1, 0) = 1;
		}
		for (size_t i = 0; i < seqTwo.length(); ++i) 
		{
			scoreMat.at(0, i + 1) = scoreMat.at(0, i) + 
									subsMat.getScore('-', seqTwo[i]);
			backtrackMat.at(0, i + 1) = 0;
		}

		//filling DP matrices
		for (size_t i = 1; i < seqOne.length() + 1; ++i)
		{
			for (size_t j = 1; j < seqTwo.length() + 1; ++j) 
			{
				float left = scoreMat.at(i, j - 1) + 
							 subsMat.getScore('-', seqTwo[j - 1]);
				float up = scoreMat.at(i - 1, j) + 
							subsMat.getScore(seqOne[i - 1], '-');
				float cross = scoreMat.at(i - 1, j - 1) + 
							  subsMat.getScore(seqOne[i - 1], seqTwo[j - 1]);

				int prev = 2;
				float score = cross;
				if (up > score)
				{
					prev = 1;
					score = up;
				}
				if (left > score)
				{
					prev = 0;
					score = left;
				}
				scoreMat.at(i, j) = score;
				backtrackMat.at(i, j) = prev;
			}
		}

		//backtrack
		int i = seqOne.length();
		int j = seqTwo.length();
		outOne.clear();
		outTwo.clear();

		while (i != 0 || j != 0) 
		{
			if(backtrackMat.at(i, j) == 1) 
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
			}
			else if (backtrackMat.at(i, j) == 0) 
			{
				outOne += '-';
				outTwo += seqTwo[j - 1];
				j -= 1;
			}
			else
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[j - 1];
				i -= 1;
				j -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
		outOne += "$";
		outTwo += "$";

		return scoreMat.at(seqOne.length(), seqTwo.length());
	}

	//Splits aligned strings into homopolymer runs (wrt to candAln)
	std::vector<std::pair<HopoMatrix::State, HopoMatrix::Observation>>
	splitBranchHopos(const std::string& candAln, const std::string& branchAln,
					 float readLikelihood)
	{
		//std::cerr << candAln << std::endl << branchAln << std::endl << std::endl;
		std::vector<std::pair<HopoMatrix::State, 
							  HopoMatrix::Observation>> result;
		size_t prevPos = 0;
		while (candAln[prevPos] == '-') ++prevPos;
		char prevNucl = candAln[prevPos];

		int runLength = 0;
		int gapLength = 0;
		int hopoCount = 0;
		for (size_t pos = prevPos + 1; pos < candAln.length(); ++pos)
		{
			if (candAln[pos] != '-') ++runLength;
			if (candAln[pos] == '-')
			{
				++gapLength;
			}
			else 
			{
				if (candAln[pos] != prevNucl)
				{
					bool leftMatch = prevPos == 0 || 
							candAln[prevPos - 1] == branchAln[prevPos - 1] ||
							branchAln[prevPos - 1] == candAln[prevPos];
					bool rightMatch = pos == candAln.length() - 1 || 
							candAln[pos] == branchAln[pos] ||
							branchAln[pos] == candAln[pos - 1];

					++hopoCount;
					//wobble
					size_t branchPrevPos = prevPos;
					if (branchAln[prevPos - 1] == candAln[prevPos]) 
						--branchPrevPos;
					size_t branchPos = pos;
					if (branchAln[pos] == candAln[pos - 1]) 
						++branchPos;

					/*if (prevPos > 0 && prevPos < branchAln.length() - 1)
					{
						std::string branchSubseq(branchAln, prevPos - 1, pos - prevPos + 2);
						std::string candSubseq(candAln, prevPos - 1, pos - prevPos + 2);
						std::cout << candSubseq << " " << std::endl 
								  << branchSubseq << std::endl 
								  << (rightMatch && leftMatch) << std::endl;
					}*/

					auto state = HopoMatrix::State(candAln, prevPos, pos);
					auto observ = HopoMatrix::strToObs(state.nucl, branchAln, 
													   branchPrevPos, branchPos);
					observ.extactMatch = leftMatch && rightMatch;
					observ.readLikelihood = readLikelihood;
					result.emplace_back(state, observ);

					//std::cout << state.length << state.nucl << " " 
					//		  << HopoMatrix::obsToStr(observ, state.nucl) << std::endl;
					/*
					if (hopoCount == 15)
					std::cerr << state.length << state.nucl << " "
							  << branchAln.substr(branchPrevPos, branchPos - branchPrevPos) << 
							  "\t" << branchAln.substr(branchPrevPos - 1, 
							  						   branchPos - branchPrevPos + 2) 
							  << std::endl;
							  */

					prevNucl = candAln[pos];
					prevPos = pos - gapLength;
					gapLength = 0;
					runLength = 0;
				}
				else
				{
					gapLength = 0;
				}
			}
		}

		return result;
	}
}

/*
float _mean(const std::vector<float>& vals)
{
	float sum = 0;
	for(float v : vals) sum += v;
	return (float)sum / vals.size();
};

float _std(const std::vector<float>& vals)
{
	float mean = _mean(vals);
	float sumStd = 0.0f;
	for (float v : vals) sumStd += (v - mean) * (v - mean);
	return sqrt(sumStd / vals.size());  
};

float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}*/


void splitBubble(const Bubble& bubble, std::vector<HopoMatrix::State>& states,
				 std::vector<HopoMatrix::ObsVector>& observations,
				 const SubstitutionMatrix& subsMatrix)
{
	for (auto& branch : bubble.branches)
	{
		std::string alnCand;
		std::string alnBranch;
		float score = pairwiseAlignment(bubble.candidate, branch, subsMatrix,
						  				alnCand, alnBranch);

		auto splitHopo = splitBranchHopos(alnCand, alnBranch, score);
		if (states.empty())
		{
			states.assign(splitHopo.size(), HopoMatrix::State());
			observations.assign(splitHopo.size(), HopoMatrix::ObsVector());
		}
		assert(states.size() == splitHopo.size());

		for (size_t i = 0; i < splitHopo.size(); ++i)
		{
			states[i] = splitHopo[i].first;
			observations[i].push_back(splitHopo[i].second);
		}
	}
}

//processes a single bubble
void HomoPolisher::polishBubble(Bubble& bubble) const
{
	//if (bubble.position != 51314) return;
	std::string prevCandidate;
	std::string curCandidate = bubble.candidate;

	std::vector<HopoMatrix::State> states;
	std::vector<HopoMatrix::ObsVector> observations;
	splitBubble(bubble, states, observations, _subsMatrix);

	const size_t MIN_HOPO = 1;
	const size_t MAX_HOPO = 10;

	typedef std::pair<double, size_t> ScorePair;
	std::string newConsensus;
	//std::cerr << bubble.position << std::endl;
	size_t candPos = 0;
	for (size_t i = 0; i < states.size(); ++i)
	{
		size_t length = states[i].length;
		if (length > 1)	//only homopolymers
		{
			/////////
			std::vector<ScorePair> scores;
			for (size_t len = std::max(MIN_HOPO, length - 1); 
				 len <= std::min(MAX_HOPO, length + 1); ++len)
			{
				/////////
				std::vector<HopoMatrix::State> newStates;
				std::vector<HopoMatrix::ObsVector> newObs;
				Bubble newBubble(bubble);
				newBubble.candidate = newBubble.candidate.substr(0, candPos) +
								  	  std::string(len, states[i].nucl) +
									  newBubble.candidate.substr(candPos + length);
				splitBubble(newBubble, newStates, newObs, _subsMatrix);

				/*
				for (auto obs : newObs[i])
				{
					std::cerr << HopoMatrix::obsToStr(obs, states[i].nucl) 
							  << "\t" << obs.extactMatch
							  << "\t" << _hopoMatrix.getObsProb(states[i], obs)
							  << std::endl;
				}
				std::cerr << std::endl;
				*/

				auto newState = HopoMatrix::State(states[i].nucl, len);
				double likelihood = this->likelihood(newState, newObs[i]);
				scores.push_back(std::make_pair(likelihood, len));
				std::cerr << likelihood << std::endl << std::endl;
			}
			std::cerr << std::endl;

			std::sort(scores.begin(), scores.end(), 
					  [](const ScorePair& p1, const ScorePair& p2)
					  {return p1.first > p2.first;});

			size_t maxRun = this->compareTopTwo(states[i].nucl, scores[0].second, 
												scores[1].second, observations[i]);
			length = scores[0].second;

			/////////
			//length = this->mostLikelyLen(states[i].nucl, observations[i]);
			std::cerr << states[i].length << states[i].nucl << " -> "
					  << length << states[i].nucl << std::endl;
		}
		candPos += states[i].length;
		newConsensus += std::string(length, states[i].nucl);

		/*
		if (length != (size_t)states[i].length)
		{
			std::cerr << bubble.position << std::endl;
			std::cerr << (int)states[i].length << states[i].nucl 
					  << " -> " << length << states[i].nucl << std::endl;
			for (auto obs : observations[i]) 
			{
				//bool outlier = obs.readLikelihood < mean - std;
				auto newState = HopoMatrix::State(states[i].nucl, length);
				std::cerr << HopoMatrix::obsToStr(obs, states[i].nucl) 
						  << "\t" << _hopoMatrix.getObsProb(states[i], obs)
						  << "\t" << _hopoMatrix.getObsProb(newState, obs) 
						  << std::endl;
			}
			std::cerr << std::endl;
		}*/
	}

	if (newConsensus != bubble.candidate)
	{
		StepInfo info;
		info.methodUsed = StepHopo;
		info.sequence = newConsensus;
		bubble.polishSteps.push_back(info);
		bubble.candidate = newConsensus;
	}
}

//likelihood of a goven state
double HomoPolisher::likelihood(HopoMatrix::State state, 
								const HopoMatrix::ObsVector& observations) const
{
	double likelihood = 0.0f;
	for (auto obs : observations)
	{
		if (obs.extactMatch)
		{
			likelihood += _hopoMatrix.getObsProb(state, obs);
		}
	}
	//likelihood += _hopoMatrix.getGenomeProb(state);
	return likelihood;
}

//for a given homopolymer and set of observations,
//computes most likely homopolymer length
size_t HomoPolisher::mostLikelyLen(char nucleotide,
								   const HopoMatrix::ObsVector& 
								   observations) const
{
	assert(!observations.empty());
	const size_t MIN_HOPO = 1;
	const size_t MAX_HOPO = 10;

	typedef std::pair<double, size_t> ScorePair;
	std::vector<ScorePair> scores;
	for (size_t len = MIN_HOPO; len <= MAX_HOPO; ++len)
	{
		auto newState = HopoMatrix::State(nucleotide, len);
		double likelihood = this->likelihood(newState, observations);
		scores.push_back(std::make_pair(likelihood, len));
		//std::cerr << likelihood << " ";
	}
	//std::cerr << std::endl;

	std::sort(scores.begin(), scores.end(), 
			  [](const ScorePair& p1, const ScorePair& p2)
			  {return p1.first > p2.first;});

	size_t maxRun = this->compareTopTwo(nucleotide, scores[0].second, 
										scores[1].second, observations);
	//return scores[0].second;
	return maxRun;
}

//Compares top two homopolimer candidates in a more precise manner
size_t HomoPolisher::compareTopTwo(char nucleotide, size_t firstChoice, 
								   size_t secondChoice,
				   				   const HopoMatrix::ObsVector& 
								   observations) const
{
	size_t choices[] = {firstChoice, secondChoice};
	HopoMatrix::ObsVector knownObs[2];

	//std::cerr << bubble.position << std::endl;
	/*
	if (firstChoice > 3 && secondChoice > 3)
	{
		std::cerr << firstChoice << nucleotide
				  << " vs " << secondChoice << nucleotide << std::endl;
		for (auto obs : observations) 
		{
			//bool outlier = obs.readLikelihood < mean - std;
			auto state1 = HopoMatrix::State(nucleotide, firstChoice);
			auto state2 = HopoMatrix::State(nucleotide, secondChoice);
			std::cerr << HopoMatrix::obsToStr(obs, nucleotide) 
					  << "\t" << _hopoMatrix.getObsProb(state1, obs)
					  << "\t" << _hopoMatrix.getObsProb(state2, obs) 
					  << std::endl;
		}
		std::cerr << std::endl;
	}*/

	for (size_t i = 0; i < 2; ++i)
	{
		auto state = HopoMatrix::State(nucleotide, choices[i]);
		knownObs[i] = _hopoMatrix.knownObservations(state);
	}
	
	//getting common known observations
	std::unordered_set<uint32_t> fstSet;
	for (auto obs : knownObs[0]) fstSet.insert(obs.id);
	std::unordered_set<uint32_t> commonSet;
	for (auto obs : knownObs[1])
	{
		if (fstSet.count(obs.id)) commonSet.insert(obs.id);
	}
	HopoMatrix::ObsVector commonObservations;
	for (auto obs : observations)
	{
		if (commonSet.count(obs.id)) commonObservations.push_back(obs);
	}

	double likelihoods[2];
	for (size_t i = 0; i < 2; ++i)
	{
		auto state = HopoMatrix::State(nucleotide, choices[i]);
		likelihoods[i] = this->likelihood(state, commonObservations);
	}

	return (likelihoods[0] > likelihoods[1]) ? choices[0] : choices[1];
}
