module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <string>
#endif

export module hardware_info;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import json;

/**
 * \namespace HardwareInfo
 * \brief Provides functionalities to retrieve and represent hardware and system information.
 *
 * The HardwareInfo namespace encapsulates functions that gather detailed information
 * about the system's hardware, compiler, operating system, and other relevant metrics.
 * It offers both human-readable string representations and structured JSON formats.
 */
export namespace HardwareInfo
{
/**
 * \brief Retrieves detailed hardware and system information as a formatted string.
 *
 * This function gathers information such as compiler details, operating system type,
 * CPU architecture, memory status, and other system-specific data. It formats this
 * information into a readable string for logging or display purposes.
 *
 * \return A string containing the formatted hardware and system information.
 */
std::string writeInfo();

/**
 * \brief Retrieves detailed hardware and system information as a JSON object.
 *
 * This function gathers the same information as `writeInfo` but structures it into
 * a JSON format, facilitating easy parsing, storage, and interoperability with other
 * systems or tools that consume JSON data.
 *
 * \return A `nlohmann::json` object containing the hardware and system information.
 */
nlohmann::json jsonInfo();
}  // namespace HardwareInfo
