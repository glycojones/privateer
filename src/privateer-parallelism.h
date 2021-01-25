
// Restraint-handling code for Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <future>
#include <mutex>
#include <queue>

#ifndef PRIVATEER_PARALLELISM_H_INCLUDED
#define PRIVATEER_PARALLELISM_H_INCLUDED


namespace privateer 
{

	namespace control
	{
        template <typename T>
        class Queue 
        {
            public:
                bool push(T const & value) 
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_queue.push(value);
                    return true;
                }

                // deletes the retrieved element, do not use for non integral types
                bool pop(T & v) 
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    if (m_queue.empty())
                        return false;
                    v = m_queue.front();
                    m_queue.pop();
                    return true;
                }

                bool empty() 
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    return m_queue.empty();
                }

                size_t queue_size(T & v)
                {
                    size_t size = m_queue.size();
                    return size; 
                }
            private:
                std::queue<T> m_queue;
                std::mutex    m_mutex;
        };
    }
	class thread_pool
    {
    public:

        thread_pool() { this->init(); }
        thread_pool(size_t nThreads) { this->init(); m_nThreads = nThreads; this->resize(nThreads); }

        // the destructor waits for all the functions in the queue to be finished
        ~thread_pool() { this->interrupt(false); }

        // get the number of running threads in the pool
        inline size_t size() const { return m_threads.size(); }

        // number of idle threads
        inline size_t n_idle() const { return ma_n_idle; }

        
        size_t n_remaining_jobs() { return this->n_jobs_in_queue(); }


        // get a specific thread
        inline std::thread & get_thread( size_t i ) { return *m_threads.at(i); }


        // restart the pool
        void restart()
        {
            this->interrupt(false); // finish all existing tasks but prevent new ones
            this->init(); // reset atomic flags
            this->resize(m_nThreads);
        }

        bool sync()
        {
            bool synced = false;
            if( (n_idle() == size()) && (n_remaining_jobs() == 0))
                return synced = true;
            else
            {
                while((n_idle() < size()) || (n_remaining_jobs() > 0))
                {
                    if( (n_idle() == size()) && (n_remaining_jobs() == 0) ) return synced = true;
                }
            }
        }

        // change the number of threads in the pool
        // should be called from one thread, otherwise be careful to not interleave, also with this->interrupt()
        void resize(size_t nThreads)
        { 
            if (!ma_kill && !ma_interrupt) 
            {
                size_t oldNThreads = m_threads.size();
                if (oldNThreads <= nThreads) {  // if the number of threads is increased

                    m_threads .resize(nThreads);
                    m_abort   .resize(nThreads);

                    for (size_t i = oldNThreads; i < nThreads; ++i)
                    {
                        m_abort[i] = std::make_shared<std::atomic<bool>>(false);
                        this->setup_thread(i);
                    }
                }
                else 
                {  // the number of threads is decreased
                    for (size_t i = oldNThreads - 1; i >= nThreads; --i)
                    {
                        *m_abort[i] = true;  // this thread will finish
                        m_threads[i]->detach();
                    }
                    {
                        // stop the detached threads that were waiting
                        std::unique_lock<std::mutex> lock(m_mutex);
                        m_cond.notify_all();
                    }
                    m_threads .resize(nThreads); // safe to delete because the threads are detached
                    m_abort   .resize(nThreads); // safe to delete because the threads have copies of shared_ptr of the flags, not originals
                }
            }
        }

        // wait for all computing threads to finish and stop all threads
        // may be called asynchronously to not pause the calling thread while waiting
        // if kill == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
        void interrupt( bool kill )
        {
            if (kill) {
                if (ma_kill) return;
                ma_kill = true;

                for (size_t i = 0, n = this->size(); i < n; ++i)
                    *m_abort[i] = true;  // command the threads to stop
            }
            else {
                if (ma_interrupt || ma_kill) return;
                ma_interrupt = true;  // give the waiting threads a command to finish
            }

            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_cond.notify_all();  // stop all waiting threads
            }
            // wait for the computing threads to finish
            for (size_t i = 0; i < m_threads.size(); ++i)
            {
                if (m_threads[i]->joinable())
                    m_threads[i]->join();
            }

            this->clear_queue();

            m_threads .clear();
            m_abort   .clear();
        }
        
        template<typename F, typename... Rest>
        auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))>
        {
            if (!ma_kill && !ma_interrupt)
            {
                auto pck = std::make_shared<std::packaged_task< decltype(f(0, rest...)) (size_t) >>(
                    std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
                );
                auto _f  = new std::function<void(size_t id)>([pck](size_t id){ (*pck)(id); });

                m_queue.push(_f);
                std::unique_lock<std::mutex> lock(m_mutex);
                m_cond.notify_one();
                return pck->get_future();
            }
            else return std::future<decltype(f(0, rest...))>();
        }


        // run the user's function that accepts argument size_t - id of the running thread. returned value is templatized
        // operator returns std::future, where the user can get the result and rethrow the catched exceptions
        template<typename F>
        auto push(F && f) ->std::future<decltype(f(0))>
        {
            if (!ma_kill && !ma_interrupt)
            {
                auto pck = std::make_shared<std::packaged_task< decltype(f(0)) (size_t) >>(std::forward<F>(f));
                auto _f  = new std::function<void(size_t id)>([pck](size_t id){ (*pck)(id); });

                m_queue.push(_f);
                std::unique_lock<std::mutex> lock(m_mutex);
                m_cond.notify_one();
                return pck->get_future();
            }
            else return std::future<decltype(f(0))>();
        }

    private:

        // deleted
        thread_pool(const thread_pool &);// = delete;
        thread_pool(thread_pool &&);// = delete;
        thread_pool & operator=(const thread_pool &);// = delete;
        thread_pool & operator=(thread_pool &&);// = delete;

        // clear all tasks
        void clear_queue()
        {
            std::function<void(size_t id)> * _f = nullptr;
            while (m_queue.pop(_f))
                delete _f; // empty the queue
        }

        size_t n_jobs_in_queue()
        {
            std::function<void(size_t id)> * _f = nullptr;
            return m_queue.queue_size(_f);
        }

        // reset all flags
        void init() { ma_n_idle = 0; ma_kill = false; ma_interrupt = false; }

        // each thread pops jobs from the queue until:
        //  - the queue is empty, then it waits (idle)
        //  - its abort flag is set (terminate without emptying the queue)
        //  - a global interrupt is set, then only idle threads terminate
        void setup_thread( size_t i )
        {
            // a copy of the shared ptr to the abort
            std::shared_ptr<std::atomic<bool>> abort_ptr(m_abort[i]);

            auto f = [this, i, abort_ptr]()
            {
                std::atomic<bool> & abort = *abort_ptr;
                std::function<void(size_t id)> * _f = nullptr;
                bool more_tasks = m_queue.pop(_f);

                while (true)
                {
                    while (more_tasks) // if there is anything in the queue
                    {
                        // at return, delete the function even if an exception occurred
                        std::unique_ptr<std::function<void(size_t id)>> func(_f);
                        (*_f)(i);
                        if (abort)
                            return; // return even if the queue is not empty yet
                        else
                            more_tasks = m_queue.pop(_f);
                    }

                    // the queue is empty here, wait for the next command
                    std::unique_lock<std::mutex> lock(m_mutex);
                    ++ma_n_idle;
                    m_cond.wait(lock, [this, &_f, &more_tasks, &abort](){ more_tasks = m_queue.pop(_f); return abort || ma_interrupt || more_tasks; });
                    --ma_n_idle;
                    if ( ! more_tasks) return; // we stopped waiting either because of interruption or abort
                }
            };

            m_threads[i].reset( new std::thread(f) );
        }
        

        // ----------  =====  ----------

        std::vector<std::unique_ptr<std::thread>>        m_threads;
        std::vector<std::shared_ptr<std::atomic<bool>>>  m_abort;
        control::Queue<std::function<void(size_t id)> *>  m_queue;

        size_t queue_size_control;

        std::atomic<bool>    ma_interrupt, ma_kill;
        std::atomic<size_t>  ma_n_idle;
        std::atomic<size_t> m_nThreads;

        std::mutex               m_mutex;
        std::condition_variable  m_cond;
    };
}

#endif
