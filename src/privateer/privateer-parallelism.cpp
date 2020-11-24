
// Restraint-handling code for Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include "privateer-parallelism.h"

template <typename T>

bool privateer::control::Queue<T>::push(T const & value) 
{
    std::unique_lock<std::mutex> lock(m_mutex);
    m_queue.push(value);
    return true;
}

// deletes the retrieved element, do not use for non integral types
template <typename T>
bool privateer::control::Queue<T>::pop(T & v) 
{
    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_queue.empty())
        return false;
    v = m_queue.front();
    m_queue.pop();
    return true;
}

template <typename T>
bool privateer::control::Queue<T>::empty() 
{
    std::unique_lock<std::mutex> lock(m_mutex);
    return m_queue.empty();
}


void privateer::thread_pool::resize(size_t nThreads)
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
void privateer::thread_pool::interrupt( bool kill )
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
auto privateer::thread_pool::push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))>
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
auto privateer::thread_pool::push(F && f) ->std::future<decltype(f(0))>
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

void privateer::thread_pool::setup_thread( size_t i )
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



void privateer::thread_pool::thread_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn, bool par){

    //internal loop
    auto int_fn = [&fn](unsigned int_start, unsigned seg_size){
        for (unsigned j = int_start; j < int_start+seg_size; j++){
            fn(j);
        }
    };

    //sequenced for
    if(!par){
        return int_fn(start, end);
    }

    //get number of available threads from already spawned thread pool.
    unsigned nb_threads_spawned = size();
    unsigned nb_threads_idle = n_idle();

    //calculate segments
    unsigned total_length = end - start;
    unsigned seg = total_length/nb_threads_idle;
    unsigned last_seg = seg + total_length%nb_threads_idle;

    //launch threads - parallel for
    auto threads_idle_vec = std::vector<std::thread>();
    auto movedindices = std::vector<unsigned>();
    threads_idle_vec.reserve(nb_threads_idle);
    for(int i = 0; i < nb_threads_spawned; i++)
    {
        unsigned current_fragment = 0;
        if(current_fragment < nb_threads_idle)
        {
            if (m_threads[i]->joinable())
            {
                unsigned current_start = seg*current_fragment;
                std::thread tempThread;
                tempThread = std::move(get_thread(i));
                tempThread = std::thread(int_fn, current_start, seg);
                threads_idle_vec.emplace_back(std::move(tempThread));
                
                movedindices.emplace_back(i);

                current_fragment++;
            }
        }
        else // sort out last fragment
        {
            if (!m_threads[i]->joinable()) // if last fragment gets a thread that is actually busy, then try next ones until we get a free thread.
            {
                continue;
            }
            else
            {
                unsigned current_start = seg*(nb_threads_idle-1);
                std::thread tempThread;
                tempThread = std::move(get_thread(i));
                tempThread = std::thread(int_fn, current_start, last_seg);
                threads_idle_vec.emplace_back(std::move(tempThread));

                movedindices.emplace_back(i);

                break; // exit the loop when we are on the last segment. 
            }
        }
    }

    for (auto& th : threads_idle_vec){
        th.join();
    }

    for(auto& index : movedindices)
    {
        setup_thread(index);
    }
}




// void privateer::thread_pool::async_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn, bool par){

//     //internal loop
//     auto int_fn = [&fn](unsigned int_start, unsigned seg_size){
//         for (unsigned j = int_start; j < int_start+seg_size; j++){
//             fn(j);
//         }
//     };
//     //sequenced for
//     if(!par){
//         return int_fn(start, end);
//     }

//     //get number of threads
//     unsigned nb_threads_hint = std::thread::hardware_concurrency();
//     unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

//     //calculate segments
//     unsigned total_length = end - start;
//     unsigned seg = total_length/nb_threads;
//     unsigned last_seg = seg + total_length%nb_threads;

//     //launch threads - parallel for
//     auto fut_vec = std::vector<std::future<void>>();
//     fut_vec.reserve(nb_threads);
//     for(int k = 0; k < nb_threads-1; ++k){
//         unsigned current_start = seg*k;
//         fut_vec.emplace_back(std::async(int_fn, current_start, seg));
//     }
//     {
//         unsigned current_start = seg*(nb_threads-1);
//         fut_vec.emplace_back(std::async(std::launch::async, int_fn, current_start, last_seg));
//     }
//     for (auto& th : fut_vec){
//         th.get();
//     }
// }
