//
// Created by anonymous authors on 2024/2/27.
//

#include "relation.h"
#include "utils.h"

#include <algorithm>



static size_t countArrayMatchesInRange(const VertexID *values, ui beginPos, ui endPos,
                                       VertexID *partMatch, const std::vector<VertexID> &toCheck) {
    if (values == nullptr || beginPos >= endPos || toCheck.empty()) return 0;
    // Candidate refinement excludes values already present in the partial match.
    if (toCheck.size() == 1) {
        VertexID v = partMatch[toCheck[0]];
        return checkExists(values, beginPos, endPos, v) != endPos ? 1 : 0;
    }

    std::vector<VertexID> checkedValues;
    checkedValues.reserve(toCheck.size());
    for (VertexID u: toCheck) checkedValues.push_back(partMatch[u]);
    std::sort(checkedValues.begin(), checkedValues.end());
    checkedValues.erase(std::unique(checkedValues.begin(), checkedValues.end()), checkedValues.end());

    size_t matched = 0;
    size_t j = 0;
    for (ui i = beginPos; i < endPos && j < checkedValues.size(); ++i) {
        while (j < checkedValues.size() && checkedValues[j] < values[i]) ++j;
        if (j < checkedValues.size() && checkedValues[j] == values[i]) {
            ++matched;
            ++j;
        }
    }
    return matched;
}

ui checkExists(const VertexID *array, const ui begin, const ui end, const ui target) {
    ui offset_begin = begin;
    ui offset_end = end;
    while (offset_end - offset_begin >= BINARY_SEARCH_THRESHOLD) {
        auto mid = (offset_begin + offset_end) / 2;
#ifndef __APPLE__
        _mm_prefetch((char *) &array[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &array[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
#endif
        if (array[mid] == target) {
            return mid;
        } else if (array[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (array[offset] == target) {
            return offset;
        }
    }

    return end;
}

size_t TrieLevel::numResults(VertexID *partMatch, const std::vector<VertexID> &toCheck,
                             const std::vector<VertexID> &largerAttrs,
                             const std::vector<VertexID> &smallerAttrs) const {
    VertexID maxCompared = 0, minCompared = (VertexID) - 1;
    for (int i = 0; i < largerAttrs.size(); ++i) {
        if (partMatch[largerAttrs[i]] > maxCompared)
            maxCompared = partMatch[largerAttrs[i]];
    }
    for (int i = 0; i < smallerAttrs.size(); ++i) {
        if (partMatch[smallerAttrs[i]] < minCompared)
            minCompared = partMatch[smallerAttrs[i]];
    }
    ui beginPos = 0, endPos = length;
    if (!largerAttrs.empty()) {
        beginPos = setBeginPos(values, length, maxCompared);
    }
    if (!smallerAttrs.empty()) {
        endPos = setEndPos(values, length, minCompared);
    }
    size_t num = endPos - beginPos;
    num -= countArrayMatchesInRange(values, beginPos, endPos, partMatch, toCheck);
    return num;
}

void TrieLevel::buildTrieFromSortedMatchesBatch(const std::vector<std::vector<VertexID>> &matches, ui start, ui end,
                                                ui depth) {
    if (matches.empty() || start >= end) {
        length = 0;
        return;
    }
    if (depth == matches[start].size() - 1) {
        length = end - start;
        delete[] values;
        values = new VertexID[length];
        for (ui idx = start; idx < end; ++idx)
            values[idx - start] = matches[idx][depth];
        return;
    }
    // Scan once and record each equal-key group.
    std::vector<std::tuple<VertexID, ui, ui>> groups;
    ui i = start;

    while (i < end) {
        VertexID currentKey = matches[i][depth];
        ui groupStart = i;

        // Find the half-open range for the current key.
        while (i < end && matches[i][depth] == currentKey) {
            ++i;
        }

        groups.emplace_back(currentKey, groupStart, i);
    }

    // Allocate one child per group.
    length = groups.size();
    delete[] values;
    values = new VertexID[length];

    // Fill the children and recursively build the next level.
    delete[] children;
    children = new TrieLevel[length];
    for (ui idx = 0; idx < length; ++idx) {
        auto [key, groupStart, groupEnd] = groups[idx];
        values[idx] = key;
        children[idx].buildTrieFromSortedMatchesBatch(matches, groupStart, groupEnd, depth + 1);
    }
}
