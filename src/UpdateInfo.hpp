#ifndef UpdateInfo_hpp
#define UpdateInfo_hpp

#include <map>
#include <string>

struct UpdateStats {
    int numTries;
    int numAccepted;
};



class UpdateInfo {

    public:
        static UpdateInfo&                  updateInfo(void)
                                                {
                                                static UpdateInfo singleUpdateInfo;
                                                return singleUpdateInfo;
                                                }
        void                                accept(void);
        void                                attempt(std::string s) { attemptedUpdateName = s; }
        void                                reject(void);
        void                                summary(void);
        
    private:
                                            UpdateInfo(void) { attemptedUpdateName = ""; }
                                            UpdateInfo(UpdateInfo& r);
        UpdateInfo&                         operator=(const UpdateInfo&);
        std::map<std::string,UpdateStats>   stats;
        std::string                         attemptedUpdateName;
};

#endif
