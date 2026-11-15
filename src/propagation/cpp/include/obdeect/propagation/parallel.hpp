#pragma once

// Row blocks of one chunk on several threads. Every block writes only its own rows and reads
// the shared inputs, so the result does not depend on the thread count (docs/propagation).

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace obdeect::propagation {

// Rows per block: large enough that a thread start is small against the work of a block.
inline constexpr std::size_t row_block = 4096;

// Runs body(begin, end) over disjoint blocks covering [0, rows) on up to `threads` threads (the
// calling thread works too). With one thread or one block, body(0, rows) runs in place. Every
// block runs even after an exception; the exception of the lowest block is rethrown after all
// threads have joined, so the error is the same for any thread count.
template <class Body> void for_row_blocks(std::size_t rows, unsigned threads, const Body &body) {
  const std::size_t blocks = (rows + row_block - 1) / row_block;
  const std::size_t workers =
      std::max<std::size_t>(1, std::min<std::size_t>(std::max(1U, threads), blocks));
  if (workers == 1) {
    body(std::size_t{0}, rows);
    return;
  }
  std::atomic<std::size_t> next{0};
  std::mutex error_mutex;
  std::size_t error_block = blocks;
  std::exception_ptr error;
  const auto work = [&] {
    for (std::size_t b = next.fetch_add(1); b < blocks; b = next.fetch_add(1)) {
      try {
        body(b * row_block, std::min(rows, (b + 1) * row_block));
      } catch (...) {
        const std::lock_guard<std::mutex> lock(error_mutex);
        if (b < error_block) {
          error_block = b;
          error = std::current_exception();
        }
      }
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(workers - 1);
  try {
    for (std::size_t w = 1; w < workers; ++w) {
      pool.emplace_back(work);
    }
  } catch (...) { // NOLINT(bugprone-empty-catch)
    // Fewer threads than asked for (e.g. the system limit): the started threads and this one
    // still cover every block.
  }
  work();
  for (auto &thread : pool) {
    thread.join();
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

} // namespace obdeect::propagation
