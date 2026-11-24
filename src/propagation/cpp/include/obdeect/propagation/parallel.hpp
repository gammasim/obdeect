#pragma once

// Row blocks of one chunk on several threads. Every block writes only its own rows and reads
// the shared inputs, so the result does not depend on the thread count
// (docs/propagation/propagation.md).
//
// RowPool keeps its worker threads for a whole run: a run calls RowPool::run several times per
// chunk, and starting threads for each call would cost more than the work of a small chunk.

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <mutex>
#include <system_error>
#include <thread>
#include <vector>

namespace obdeect::propagation {

// Rows per block: large enough that a block's work is large against the synchronisation.
inline constexpr std::size_t row_block = 4096;

// Number of blocks covering [0, rows).
[[nodiscard]] inline std::size_t row_block_count(std::size_t rows) noexcept {
  return (rows + row_block - 1) / row_block;
}

// A fixed set of worker threads (threads - 1 of them; the calling thread works too). run()
// covers [0, rows) with row blocks. With one thread or one block, body(0, rows) runs in place.
// Every block runs even after an exception; the exception of the lowest block is rethrown after
// all workers have finished the run, so the error does not depend on the thread count. A
// RowPool is used by one thread at a time.
class RowPool {
public:
  // threads: at least 1. Fewer workers are started if the system refuses a thread.
  explicit RowPool(unsigned threads) {
    const unsigned wanted = std::max(1U, threads);
    workers_.reserve(wanted - 1U);
    try {
      for (unsigned i = 1; i < wanted; ++i) {
        workers_.emplace_back([this] { worker_loop(); });
      }
    } catch (const std::system_error &) { // NOLINT(bugprone-empty-catch)
      // Fewer threads than asked for: the started workers and the caller still cover every block.
    }
  }

  ~RowPool() {
    {
      const std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    work_ready_.notify_all();
    for (auto &worker : workers_) {
      worker.join();
    }
  }

  RowPool(const RowPool &) = delete;
  RowPool &operator=(const RowPool &) = delete;
  RowPool(RowPool &&) = delete;
  RowPool &operator=(RowPool &&) = delete;

  // Number of threads that run a block (the workers and the caller).
  [[nodiscard]] unsigned threads() const noexcept {
    return static_cast<unsigned>(workers_.size()) + 1U;
  }

  // Runs body(begin, end) over disjoint blocks covering [0, rows).
  template <class Body> void run(std::size_t rows, const Body &body) {
    const std::size_t blocks = row_block_count(rows);
    if (workers_.empty() || blocks <= 1) {
      body(std::size_t{0}, rows);
      return;
    }
    const auto thunk = [](const void *context, std::size_t begin, std::size_t end) {
      (*static_cast<const Body *>(context))(begin, end);
    };
    dispatch(thunk, &body, rows, blocks);
  }

private:
  using Thunk = void (*)(const void *, std::size_t, std::size_t);

  struct Job {
    Thunk thunk{};
    const void *context{};
    std::size_t rows{};
    std::size_t blocks{};
  };

  void dispatch(Thunk thunk, const void *context, std::size_t rows, std::size_t blocks) {
    {
      const std::lock_guard<std::mutex> lock(mutex_);
      job_ = Job{thunk, context, rows, blocks};
      next_block_ = 0;
      error_block_ = blocks;
      error_ = nullptr;
      finished_ = 0;
      ++generation_;
    }
    work_ready_.notify_all();
    work();
    std::exception_ptr error;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      work_done_.wait(lock, [this] { return finished_ == workers_.size(); });
      error = error_;
    }
    if (error) {
      std::rethrow_exception(error);
    }
  }

  // Claims blocks of the current job until none is left. Runs on the caller and on the workers.
  // job_ was written under mutex_ before the workers were woken, so reading it here is ordered.
  void work() {
    const Job job = job_;
    for (std::size_t b = next_block_.fetch_add(1); b < job.blocks; b = next_block_.fetch_add(1)) {
      try {
        job.thunk(job.context, b * row_block, std::min(job.rows, (b + 1) * row_block));
      } catch (...) {
        const std::lock_guard<std::mutex> lock(mutex_);
        if (b < error_block_) {
          error_block_ = b;
          error_ = std::current_exception();
        }
      }
    }
  }

  void worker_loop() {
    std::size_t seen = 0;
    while (true) {
      {
        std::unique_lock<std::mutex> lock(mutex_);
        work_ready_.wait(lock, [this, seen] { return stop_ || generation_ != seen; });
        if (stop_) {
          return;
        }
        seen = generation_;
      }
      work();
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        ++finished_;
      }
      work_done_.notify_one();
    }
  }

  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable work_ready_;
  std::condition_variable work_done_;
  bool stop_{false};
  std::size_t generation_{0};
  std::size_t finished_{0};
  Job job_{};
  std::atomic<std::size_t> next_block_{0};
  std::size_t error_block_{0};
  std::exception_ptr error_;
};

// Runs body(begin, end) over disjoint blocks covering [0, rows) on up to `threads` threads, with
// a pool of its own for this call. A run with several calls per chunk keeps one RowPool instead.
template <class Body> void for_row_blocks(std::size_t rows, unsigned threads, const Body &body) {
  const std::size_t blocks = row_block_count(rows);
  const auto workers = static_cast<unsigned>(
      std::max<std::size_t>(1, std::min<std::size_t>(std::max(1U, threads), blocks)));
  RowPool pool(workers);
  pool.run(rows, body);
}

} // namespace obdeect::propagation
