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

#ifndef ATOMICALLYGROWINGVECTOR_H
#define ATOMICALLYGROWINGVECTOR_H

#include <vector>
#include <atomic>
#include <stdexcept>
#include <algorithm>

/**
 * Atomically growing vector compatible with std::vector.
 *
 * Mostly simulating C++11 behavior.
 * Additional back_insert method allows appending multiple consecutive items atomically.
 * It behaves similar to insert(end(), ...).
 *
 * Concurrent calls to back_insert, push_back and emplace_back are valid.
 * All other modifying accesses need to be synchronized externally.
 * To preserve validity of the returned references and iterators, automatic reallocation is disabled.
 * Therefore, reserve needs to be called in advance to provide sufficient storage for the given use-case.
 * Trying to insert elements once the pre-allocated storage is exhausted will throw std::length_error and will not modify the object's state.
 *
 * Implemented using std::vector and one additional std::atomic<std::vector<T>::iterator>.
 */
template<class T, class Allocator=std::allocator<T>>
class AtomicallyGrowingVector
{
public:
    typedef std::vector<T, Allocator> VectorType;

    //member types
    using value_type             = typename VectorType::value_type;
    using allocator_type         = typename VectorType::allocator_type;
    using size_type              = typename VectorType::size_type;
    using difference_type        = typename VectorType::difference_type;
    using reference              = typename VectorType::reference;
    using const_reference        = typename VectorType::const_reference;
    using pointer                = typename VectorType::pointer;
    using const_pointer          = typename VectorType::const_pointer;
    using iterator               = typename VectorType::iterator;
    using const_iterator         = typename VectorType::const_iterator;
    using reverse_iterator       = typename VectorType::reverse_iterator;
    using const_reverse_iterator = typename VectorType::const_reverse_iterator;

private:
    std::vector<T, Allocator> internalVector;
    std::atomic<iterator> atomicEnd;

    /// throws std::length_error if not sufficiently pre-allocated
    iterator grow(size_t numElements=1) {
        iterator oldEnd = atomicEnd.load(std::memory_order_acquire);

        if (numElements == 0)
            return oldEnd;

        iterator newEnd;
        do
        {
            newEnd = oldEnd+numElements;
            if (newEnd > internalVector.end())
                throw std::length_error("AtomicallyGrowingVector::grow(): not enough pre-allocated storage available.");
        }
        while (!atomicEnd.compare_exchange_weak(oldEnd, newEnd, std::memory_order_acq_rel, std::memory_order_acquire));
        return oldEnd;
    }

    /// shrinking below 0 elements is undefined behavior
    void shrink(size_t numElements=1) {
        if (numElements == 0)
            return;

        iterator oldEnd = atomicEnd.load(std::memory_order_acquire);
        while (!atomicEnd.compare_exchange_weak(oldEnd, oldEnd-numElements, std::memory_order_acq_rel, std::memory_order_acquire));
    }

public:
    template<typename... Ts> AtomicallyGrowingVector(Ts... parameters) : internalVector{parameters...}, atomicEnd{internalVector.end()} {}

    AtomicallyGrowingVector(const AtomicallyGrowingVector& other) : internalVector{other.internalVector}, atomicEnd{begin()+other.size()} {}
    AtomicallyGrowingVector(AtomicallyGrowingVector&& other) : internalVector{std::move(other.internalVector)}, atomicEnd{other.atomicEnd.load(std::memory_order_acquire)} {}

    //assignment
    AtomicallyGrowingVector& operator=(const AtomicallyGrowingVector& other) {
        if (this == &other)
            return *this;

        internalVector = other.internalVector;
        atomicEnd.store(begin()+other.size(), std::memory_order_release);
        return *this;
    }
    AtomicallyGrowingVector& operator=(AtomicallyGrowingVector&& other) {
        if (this == &other)
            return *this;

        internalVector = std::move(other.internalVector);
        atomicEnd.store(other.atomicEnd.load(std::memory_order_acquire), std::memory_order_release);
        return *this;
    }
    AtomicallyGrowingVector& operator=(std::initializer_list<T> ilist) {
        *this = AtomicallyGrowingVector{ilist};
        return *this;
    }

    void assign(size_type count, const T& value) {
        internalVector.assign(count, value);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }
    template<class InputIt>
    void assign(InputIt first, InputIt last) {
        internalVector.assign(first, last);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }
    void assign(std::initializer_list<T> ilist) {
        internalVector.assign(ilist);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }

    allocator_type get_allocator() const {
        return internalVector.get_allocator();
    }

    //element access
    reference at(size_type pos) {
        reference item = internalVector.at(pos);
        if (item >= atomicEnd.load(std::memory_order_acquire))
            throw std::out_of_range("AtomicallyGrowingVector::at(): requested element is out of range.");
        return item;
    }
    const_reference at(size_type pos) const {
        const_reference item = internalVector.at(pos);
        if (item >= atomicEnd.load(std::memory_order_acquire))
            throw std::out_of_range("AtomicallyGrowingVector::at(): requested element is out of range.");
        return item;
    }

    reference operator[](size_type pos) {
        return internalVector.operator[](pos);
    }
    const_reference operator[](size_type pos) const {
        return internalVector.operator[](pos);
    }

    reference front() {
        return internalVector.front();
    }
    const_reference front() const {
        return internalVector.front();
    }

    reference back() {
        return --atomicEnd.load(std::memory_order_acquire);
    }
    const_reference back() const {
        return --atomicEnd.load(std::memory_order_acquire);
    }

    T* data() noexcept {
        return internalVector.data();
    }
    const T* data() const noexcept {
        return internalVector.data();
    }

    //iterators
    iterator begin() noexcept {
        return internalVector.begin();
    }
    const_iterator begin() const noexcept {
        return internalVector.begin();
    }
    const_iterator cbegin() const noexcept {
        return internalVector.cbegin();
    }

    iterator end() noexcept {
        return atomicEnd.load(std::memory_order_acquire);
    }
    const_iterator end() const noexcept {
        return atomicEnd.load(std::memory_order_acquire);
    }
    const_iterator cend() const noexcept {
        return atomicEnd.load(std::memory_order_acquire);
    }

    reverse_iterator rbegin() noexcept {
        return internalVector.rbegin();
    }
    const_reverse_iterator rbegin() const noexcept {
        return internalVector.rbegin();
    }
    const_reverse_iterator crbegin() const noexcept {
        return internalVector.crbegin();
    }

    reverse_iterator rend() noexcept {
        return back();
    }
    const_reverse_iterator rend() const noexcept {
        return back();
    }
    const_reverse_iterator crend() const noexcept {
        return back();
    }

    //capacity
    bool empty() const noexcept {
        return begin() == end();
    }

    size_type size() const noexcept {
        return std::distance(begin(), end());
    }

    size_type max_size() const noexcept {
        return internalVector.max_size();
    }

    void reserve(size_type new_cap) {
        const size_type pos = size();
        internalVector.reserve(new_cap);
        internalVector.resize(internalVector.capacity());
        atomicEnd.store(begin()+pos, std::memory_order_release);
    }

    size_type capacity() const noexcept {
        //the internal vector's size is the capacity of the atomically growing vector
        return internalVector.size();
    }

    void shrink_to_fit() {
        const size_type pos = size();
        internalVector.shrink_to_fit();
        atomicEnd.store(begin()+pos, std::memory_order_release);
    }

    //modifiers
    void clear() noexcept {
        internalVector.clear();
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }

    iterator insert(const_iterator pos, const T& value) {
        const iterator oldEnd = grow();
        std::move_backward(pos, oldEnd, oldEnd+1);
        *pos = value;
        return pos;
    }
    iterator insert(const_iterator pos, T&& value) {
        const iterator oldEnd = grow();
        std::move_backward(pos, oldEnd, oldEnd+1);
        *pos = std::move(value);
        return pos;
    }
    iterator insert(const_iterator pos, size_type count, const T& value) {
        const iterator oldEnd = grow(count);
        std::move_backward(pos, oldEnd, oldEnd+count);
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos)
            *insertPos = value;
        return pos;
    }
    template<class InputIt>
    iterator insert(const_iterator pos, InputIt first, InputIt last) {
        size_type count = std::distance(first, last);
        const iterator oldEnd = grow(count);
        std::move_backward(pos, oldEnd, oldEnd+count);
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos, ++first)
            *insertPos = *first;
        return pos;
    }
    iterator insert(const_iterator pos, std::initializer_list<T> ilist) {
        size_type count = ilist.size();
        const iterator oldEnd = grow(count);
        std::move_backward(pos, oldEnd, oldEnd+count);
        typename std::initializer_list<T>::const_iterator first = ilist.begin();
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos, ++first)
            *insertPos = *first;
        return pos;
    }

    //atomic insertion at the end of the vector
    iterator back_insert(const T& value) {
        const iterator pos = grow();
        *pos = value;
        return pos;
    }
    iterator back_insert(T&& value) {
        const iterator pos = grow();
        *pos = std::move(value);
        return pos;
    }
    iterator back_insert(size_type count, const T& value) {
        const iterator pos = grow(count);
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos)
            *insertPos = value;
        return pos;
    }
    template<class InputIt>
    iterator back_insert(InputIt first, InputIt last) {
        size_type count = std::distance(first, last);
        const iterator pos = grow(count);
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos, ++first)
            *insertPos = *first;
        return pos;
    }
    iterator back_insert(std::initializer_list<T> ilist) {
        size_type count = ilist.size();
        const iterator pos = grow(count);
        typename std::initializer_list<T>::const_iterator ilistIter = ilist.begin();
        for (iterator insertPos=pos; insertPos!=pos+count; ++insertPos, ++ilistIter)
            *insertPos = *ilistIter;
        return pos;
    }

    template<class... Args>
    iterator emplace(const_iterator pos, Args&&... args) {
        const iterator tmp = grow();
        std::move_backward(pos, tmp, tmp+1);
        pos->~T();
        internalVector.get_allocator().construct(&*pos, std::forward<Args>(args)...);
        return pos;
    }

    iterator erase(const_iterator pos) {
        shrink();
        return internalVector.erase(pos);
    }

    iterator erase(const_iterator first, const_iterator last) {
        shrink(std::distance(first, last));
        return internalVector.erase(first, last);
    }

    void push_back(const T& value) {
        *grow() = value;
    }
    void push_back(T&& value) {
        *grow() = std::move(value);
    }

    template<class... Args>
    void emplace_back(Args&&... args) {
        const iterator tmp = grow();
        tmp->~T();
        internalVector.get_allocator().construct(&*tmp, std::forward<Args>(args)...);
    }

    void pop_back() {
        shrink();
    }

    void resize(size_type count) {
        internalVector.resize(count);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }
    void resize(size_type count, const value_type& value) {
        internalVector.resize(count, value);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }

    void swap(AtomicallyGrowingVector& other) {
        internalVector.swap(other);
        iterator tmp = atomicEnd.exchange(other.atomicEnd, std::memory_order_acq_rel);
        other.atomicEnd.store(tmp, std::memory_order_release);
    }
    void swap(VectorType& other) {
        internalVector.resize(size());
        internalVector.swap(other);
        atomicEnd.store(internalVector.end(), std::memory_order_release);
    }
};

#endif // ATOMICALLYGROWINGVECTOR_H
