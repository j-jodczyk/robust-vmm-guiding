/*
    This file is part of the implementation of the SIGGRAPH 2020 paper
    "Robust Fitting of Parallax-Aware Mixtures for Path Guiding",
    as well as the updated implementation of the ACM TOG 2019 paper
    "Volume Path Guiding Based on Zero-Variance Random Walk Theory".
    The implementation extends Mitsuba, a physically based rendering system.

    Copyright (c) 2020 Lukas Ruppert, Sebastian Herholz.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RANGE_H
#define RANGE_H

#include <mitsuba/mitsuba.h>

#include <mitsuba/core/fstream.h>

#include <string>
#include <ostream>
#include <fstream>
#include <vector>
#include <type_traits>

GUIDING_NAMESPACE_BEGIN

/**
 * convenience class for range-based for loops over subsets of containers
 */
template <typename TContainer>
class Range {
public:
    typedef TContainer container_type;
    typedef typename TContainer::value_type value_type;
    typedef typename TContainer::iterator iterator;
private:
    iterator m_begin;
    iterator m_end;

public:
    FINLINE iterator begin() const {
        return m_begin;
    }
    FINLINE iterator end() const {
        return m_end;
    }

    Range() = default;

    Range(TContainer& container) : Range{container.begin(), container.end()} {}
    Range(const iterator& begin, const iterator& end) : m_begin{begin}, m_end{end} {}

    Range(TContainer& container, typename TContainer::difference_type offsetMin, typename TContainer::difference_type count)
        : Range{container.begin(), offsetMin, count} {}
    Range(const iterator& begin, typename TContainer::difference_type offsetMin, typename TContainer::difference_type count)
        : m_begin{begin+offsetMin}, m_end{m_begin+count} {}

    FINLINE value_type& operator[](ptrdiff_t offset) {
        return *(m_begin+offset);
    }

    FINLINE const value_type& operator[](ptrdiff_t offset) const {
        return *(m_begin+offset);
    }

    FINLINE size_t size() const {
        return static_cast<size_t>(std::distance(m_begin, m_end));
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Range[size: " << size() << "]\n";
        return oss.str();
    }

    template<typename T=value_type, typename std::enable_if<std::is_base_of<SerializableObject, T>::value>::type* = 0>
    void serialize(const std::string& filename) const {
        //serializable using Mitsuba's Stream
        ref<FileStream> fileStream = new FileStream(filename, FileStream::ETruncWrite);
        const size_t numElements = size();
        fileStream->writeSize(numElements);

        auto it = begin();
        for (size_t i=0; i<numElements; ++i, ++it)
            it->serialize(fileStream);

        fileStream->close();
    }

    template<typename T=value_type, typename std::enable_if<std::is_trivially_copyable<T>::value>::type* = 0>
    void serialize(const std::string& filename) const {
        //serializable as plain data
        std::ofstream file;
        file.open(filename, std::ios::binary);
        const size_t numElements = size();
        file.write(reinterpret_cast<const char*>(&numElements), sizeof(numElements));
        auto it = begin();
        for (size_t i=0; i<numElements; ++i, ++it)
            file.write(reinterpret_cast<const char*>(&*it), sizeof(value_type));
        file.close();
    }

    template<typename T=value_type, typename std::enable_if<std::is_base_of<SerializableObject, T>::value>::type* = 0>
    static std::vector<value_type> deserialize(const std::string& filename) {
        //serializable using Mitsuba's Stream
        ref<FileStream> fileStream = new FileStream(filename, FileStream::EReadOnly);
        size_t numElements = fileStream->readSize();

        TContainer container{numElements};
        for (size_t i=0; i<numElements; ++i)
            container.emplace_back(fileStream);

        fileStream->close();
        return container;
    }

    template<typename T=value_type, typename std::enable_if<std::is_trivially_copyable<T>::value>::type* = 0>
    static TContainer deserialize(const std::string& filename) {
        //serializable as plain data
        std::ifstream file;
        file.open(filename, std::ios::binary);
        size_t numElements;
        file.read(reinterpret_cast<char*>(&numElements), sizeof(numElements));
        TContainer container{numElements};
        file.read(reinterpret_cast<char*>(container.data()), sizeof(value_type)*numElements);
        file.close();
        return container;
    }
};

///to simplify the definition of ranges for arbitrary containers, including other ranges
template<typename TContainer>
class Range<Range<TContainer>> : public Range<TContainer>
{
    template<typename... Ts> Range(Ts... parameters) : Range<TContainer>(parameters...) {}
};

GUIDING_NAMESPACE_END

#endif // RANGE_H
