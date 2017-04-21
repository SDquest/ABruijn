//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <limits>

struct FastaRecord
{
	class Id
	{
	public:
		explicit Id(uint32_t id): _id(id) {}

		bool operator==(const Id& other) const
			{return _id == other._id;}
		bool operator!=(const Id& other) const
			{return !(*this == other);}
		bool operator<(const Id& other) const
			{return _id < other._id;}

		Id rc() const		//reverse complement 
			{return Id(_id + 1 - (_id % 2) * 2);}
		bool strand() const		//true = positive, false = negative
			{return !(_id % 2);}
		size_t hash() const 
			{return 0x9ddfea08eb382d69ULL * (size_t)_id;}
		int signedId() const
			{return (_id % 2) ? -((int)_id + 1) / 2 : (int)_id / 2 + 1;}

		friend std::ostream& operator << (std::ostream& stream, const Id& id)
		{
			stream << std::to_string(id._id);
			return stream;
		}
		
		friend std::istream& operator >> (std::istream& stream, Id& id)
		{
			std::string buffer;
			stream >> buffer;
			id._id = std::stoi(buffer);
			return stream;
		}

	private:
		uint32_t _id;
	};
	static const Id ID_NONE; 
	typedef std::tuple<Id, Id> IdPair;

	class DnaRepr
	{
	public:
		explicit DnaRepr() {}
		explicit DnaRepr(const std::string& string):
			_representation(string.begin(), string.end()) {}

		size_t length() const {return _representation.size();}
		char at(size_t index) const {return _representation[index];}
		DnaRepr substr(size_t start, size_t length) const 
		{
			if (start + length > _representation.size())
			{
				length = _representation.size() - start;
			}

			DnaRepr substr;
			std::copy(_representation.begin() + start,
					  _representation.begin() + start + length,
					  std::back_inserter(substr._representation));
			return substr;
		}
		void push_back(char c) {_representation.push_back(c);}
		std::string str() 
		{
			return std::string(_representation.begin(), _representation.end());
		}
		void swap(DnaRepr& other) {_representation.swap(other._representation);}

	private:
		std::vector<char> _representation;
	};

	FastaRecord(): id(ID_NONE) {}
	FastaRecord(const DnaRepr& sequence, const std::string& description,
				Id id):
		id(id), sequence(sequence), description(description)
	{
	}

	FastaRecord(const FastaRecord& other):
		id(other.id), sequence(other.sequence), 
		description(other.description) {}

	FastaRecord(FastaRecord&& other):
		id (other.id)
	{
		*this = std::move(other);
	}

	FastaRecord& operator=(const FastaRecord& other)
	{
		id = other.id;
		sequence = other.sequence;
		description = other.description;
		return *this;
	}

	FastaRecord& operator=(FastaRecord&& other)
	{
		id = other.id;
		sequence.swap(other.sequence);
		description.swap(other.description);
		return *this;
	}
	
	Id id;
	DnaRepr sequence;
	std::string description;		
};

namespace std
{
	template <>
	struct hash<FastaRecord::Id> 
	{
		size_t operator() (const FastaRecord::Id& h) const throw() 
		{
			 return h.hash();
		}
	};

	template <>
	struct hash<FastaRecord::IdPair> 
	{
		 size_t operator()(const FastaRecord::IdPair& k) const
		 {
			size_t lhs = std::get<0>(k).hash();
			size_t rhs = std::get<1>(k).hash();
			lhs ^= rhs + 0x9ddfea08eb382d69ULL + (lhs << 6) + (lhs >> 2);
			return lhs;
		 }
	};
}


class SequenceContainer
{
public:
	typedef std::unordered_map<FastaRecord::Id, 
							   FastaRecord> SequenceIndex;

	SequenceContainer() {}

	static void writeFasta(const std::vector<FastaRecord>& records,
						   const std::string& fileName);

	const SequenceIndex& getIndex() const
	{
		return _seqIndex;
	}

	const FastaRecord::DnaRepr& getSeq(FastaRecord::Id readId) const
	{
		return _seqIndex.at(readId).sequence;
	}

	int32_t seqLen(FastaRecord::Id readId) const
	{
		return _seqIndex.at(readId).sequence.length();
	}
	std::string seqName(FastaRecord::Id readId) const
	{
		return _seqIndex.at(readId).description;
	}
	void readFasta(const std::string& filename);
	const FastaRecord&  addSequence(const FastaRecord::DnaRepr& sequence, 
									const std::string& description);

private:

	size_t 	getSequences(std::vector<FastaRecord>& record, 
						 const std::string& fileName);
	size_t 	getSequencesWithComplements(std::vector<FastaRecord>& record, 
										const std::string& fileName);
	void 	validateSequence(const std::string& sequence);
	void 	validateHeader(std::string& header);

	SequenceIndex _seqIndex;
	static size_t g_nextSeqId;
};
