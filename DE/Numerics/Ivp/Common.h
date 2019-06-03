#pragma once

#include <utility>

#include "../Arg/ArgN.h"
#include "../Utils/Utils.h"
#include "../Common/Common.h"


namespace Num
{
    namespace Ivp
    {
		template<class Argument, class Value>
        using Node = std::pair<Argument, Value>;
    }
}