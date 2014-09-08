(defproject me.lizier/jidt "@VERSION@"
  :description "Java Information Dynamics Toolkit (JIDT)"
  :url "https://code.google.com/p/information-dynamics-toolkit/"
  :mailing-list {:name "jidt-discuss"
		 :archive "http://groups.google.com/group/jidt-discuss"
		 :post "jidt-discuss@googlegroups.com"
		 :subscribe "jidt-discuss+subscribe@googlegroups.com."
		 :unsubscribe "jidt-discuss+unsubscribe@googlegroups.com"}
  :license
    {
      :name "GNU GPL v3"
      :url "http://www.gnu.org/licenses/gpl.html"
      :distribution :repo
    }
  :dependencies []
  :java-source-paths
    ["../../../java/source"]
  :javac-options ["-target" "1.6" "-source" "1.6" "-Xlint:-options"]
  :aot :all
  :omit-source true
  :signing {:gpg-key "joseph.lizier@gmail.com"}
)
