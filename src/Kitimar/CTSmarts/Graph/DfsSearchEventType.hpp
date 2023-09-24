#pragma once

namespace Kitimar::CTSmarts {

    enum class DfsSearchEventType
    {
        Invalid,
        VisitVertex, // index = vertex index, flag = isNewComponent
        VisitEdge, // index = edge index, flag = isClosure
        BacktrackVertex, // index = vertex index, flag = isEndComponent
        BacktrackEdge // index = edge index, flag = isClosure
    };

}
