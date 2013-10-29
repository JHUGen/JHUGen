


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>usercode/snowmass2013/src/RooSpinZero_5D_ZH.cc at master · awhitbeck/usercode · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="github-fe137-cp1-prd.iad.github.net">
    <meta name="ruby" content="ruby 1.9.3p194-tcs-github-tcmalloc (0e75de19f8) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035/">
    <link rel="xhr-socket" href="/_sockets" />
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="83E13811:27E9:7B1EBC:526FCB6F" name="octolytics-dimension-request_id" />
    

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="Gv4ig7I1fW1DzXOM63hxLUEyOcjxStxMG2YZPviDhdE=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-32b82bb8400eaead523f5ac84d94f34d33524ef7.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-10157052ff99006b9ff1a03c3407021933da16a0.css" media="all" rel="stylesheet" type="text/css" />
    

    

      <script src="https://github.global.ssl.fastly.net/assets/frameworks-5da1cc8ea24efcf7bd73fdfdaef699713c11f597.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-d67069913b16ce39135c781a96daa58d43bd6606.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="743fbe45095908ba654d41e39c418ec1">

        <link data-pjax-transient rel='permalink' href='/awhitbeck/usercode/blob/d2f33cc63262484aea017ce3bdbb4e624f2b38b2/snowmass2013/src/RooSpinZero_5D_ZH.cc'>
  <meta property="og:title" content="usercode"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/awhitbeck/usercode"/>
  <meta property="og:image" content="https://github.global.ssl.fastly.net/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="usercode - CMS UserCode/Whitbeck -- imported from CVS (http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Whitbeck/)"/>

  <meta name="description" content="usercode - CMS UserCode/Whitbeck -- imported from CVS (http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Whitbeck/)" />

  <meta content="4173436" name="octolytics-dimension-user_id" /><meta content="awhitbeck" name="octolytics-dimension-user_login" /><meta content="11638078" name="octolytics-dimension-repository_id" /><meta content="awhitbeck/usercode" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="11638078" name="octolytics-dimension-repository_network_root_id" /><meta content="awhitbeck/usercode" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/awhitbeck/usercode/commits/master.atom" rel="alternate" title="Recent Commits to usercode:master" type="application/atom+xml" />

  </head>


  <body class="logged_out  env-production  vis-public  page-blob">
    <div class="wrapper">
      
      
      


      
      <div class="header header-logged-out">
  <div class="container clearfix">

    <a class="header-logo-wordmark" href="https://github.com/">
      <span class="mega-octicon octicon-logo-github"></span>
    </a>

    <div class="header-actions">
        <a class="button primary" href="/join">Sign up</a>
      <a class="button signin" href="/login?return_to=%2Fawhitbeck%2Fusercode%2Fblob%2Fmaster%2Fsnowmass2013%2Fsrc%2FRooSpinZero_5D_ZH.cc">Sign in</a>
    </div>

    <div class="command-bar js-command-bar  in-repository">

      <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
        <li class="features"><a href="/features">Features</a></li>
          <li class="enterprise"><a href="https://enterprise.github.com/">Enterprise</a></li>
          <li class="blog"><a href="/blog">Blog</a></li>
      </ul>
        <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    
      data-repo="awhitbeck/usercode"
      data-branch="master"
      data-sha="dba4ca19dd511a2ee9acdaa1c0ed66e4ac29cf10"
  >

    <input type="hidden" name="nwo" value="awhitbeck/usercode" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
    </div>

  </div>
</div>


      


          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">


  <li>
  <a href="/login?return_to=%2Fawhitbeck%2Fusercode"
  class="minibutton with-count js-toggler-target star-button entice tooltipped upwards"
  title="You must be signed in to use this feature" rel="nofollow">
  <span class="octicon octicon-star"></span>Star
</a>
<a class="social-count js-social-count" href="/awhitbeck/usercode/stargazers">
  0
</a>

  </li>

    <li>
      <a href="/login?return_to=%2Fawhitbeck%2Fusercode"
        class="minibutton with-count js-toggler-target fork-button entice tooltipped upwards"
        title="You must be signed in to fork a repository" rel="nofollow">
        <span class="octicon octicon-git-branch"></span>Fork
      </a>
      <a href="/awhitbeck/usercode/network" class="social-count">
        0
      </a>
    </li>
</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author">
            <a href="/awhitbeck" class="url fn" itemprop="url" rel="author"><span itemprop="title">awhitbeck</span></a>
          </span>
          <span class="repohead-name-divider">/</span>
          <strong><a href="/awhitbeck/usercode" class="js-current-repository js-repo-home-link">usercode</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container ">

        <div class="repository-sidebar">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped leftwards" title="Code">
        <a href="/awhitbeck/usercode" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /awhitbeck/usercode">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped leftwards" title="Issues">
          <a href="/awhitbeck/usercode/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-gotokey="i" data-selected-links="repo_issues /awhitbeck/usercode/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped leftwards" title="Pull Requests"><a href="/awhitbeck/usercode/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /awhitbeck/usercode/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/awhitbeck/usercode/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_wiki /awhitbeck/usercode/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/awhitbeck/usercode/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /awhitbeck/usercode/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/awhitbeck/usercode/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /awhitbeck/usercode/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/awhitbeck/usercode/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /awhitbeck/usercode/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


  </div>
</div>

            <div class="only-with-full-nav">
              

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/awhitbeck/usercode.git" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/awhitbeck/usercode.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/awhitbeck/usercode" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/awhitbeck/usercode" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <span class="octicon help tooltipped upwards" title="Get help on which URL is right for you.">
    <a href="https://help.github.com/articles/which-remote-url-should-i-use">
    <span class="octicon octicon-question"></span>
    </a>
  </span>
</p>



              <a href="/awhitbeck/usercode/archive/master.zip"
                 class="minibutton sidebar-button"
                 title="Download this repository as a zip file"
                 rel="nofollow">
                <span class="octicon octicon-cloud-download"></span>
                Download ZIP
              </a>
            </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:689d9fcc21d61def566a6f0db28ecef0 -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/awhitbeck/usercode/find/master" data-pjax data-hotkey="t" class="js-show-file-finder" style="display:none">Show File Finder</a>

<div class="file-navigation">
  
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    role="button" aria-label="Switch branches or tags" tabindex="0">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/blob/master/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/whitbeck_v1/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="whitbeck_v1"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="whitbeck_v1">whitbeck_v1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/superMELAresultsHCP/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="superMELAresultsHCP"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="superMELAresultsHCP">superMELAresultsHCP</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/quasi_stable_V00/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="quasi_stable_V00"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="quasi_stable_V00">quasi_stable_V00</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/preLHCacceptance/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="preLHCacceptance"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="preLHCacceptance">preLHCacceptance</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/ilc_angular_v1/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="ilc_angular_v1"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="ilc_angular_v1">ilc_angular_v1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/awhitbeck/usercode/tree/HCP_spin_unblinding/snowmass2013/src/RooSpinZero_5D_ZH.cc"
                 data-name="HCP_spin_unblinding"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="HCP_spin_unblinding">HCP_spin_unblinding</a>
            </div> <!-- /.select-menu-item -->
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/awhitbeck/usercode" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">usercode</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/awhitbeck/usercode/tree/master/snowmass2013" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">snowmass2013</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/awhitbeck/usercode/tree/master/snowmass2013/src" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">src</span></a></span><span class="separator"> / </span><strong class="final-path">RooSpinZero_5D_ZH.cc</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="snowmass2013/src/RooSpinZero_5D_ZH.cc" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  <div class="commit commit-loader file-history-tease js-deferred-content" data-url="/awhitbeck/usercode/contributors/master/snowmass2013/src/RooSpinZero_5D_ZH.cc">
    Fetching contributors…

    <div class="participation">
      <p class="loader-loading"><img alt="Octocat-spinner-32-eaf2f5" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" /></p>
      <p class="loader-error">Cannot retrieve contributors at this time</p>
    </div>
  </div>

<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">executable file</span>
          <span>533 lines (405 sloc)</span>
        <span>14.688 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
              <a class="minibutton disabled js-entice" href=""
                 data-entice="You must be signed in to make or propose changes">Edit</a>
          <a href="/awhitbeck/usercode/raw/master/snowmass2013/src/RooSpinZero_5D_ZH.cc" class="button minibutton " id="raw-url">Raw</a>
            <a href="/awhitbeck/usercode/blame/master/snowmass2013/src/RooSpinZero_5D_ZH.cc" class="button minibutton ">Blame</a>
          <a href="/awhitbeck/usercode/commits/master/snowmass2013/src/RooSpinZero_5D_ZH.cc" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
          <a class="minibutton danger empty-icon js-entice" href=""
             data-entice="You must be signed in and on a branch to make or propose changes">
          Delete
        </a>
      </div><!-- /.actions -->

    </div>
        <div class="blob-wrapper data type-c js-blob-data">
        <table class="file-code file-diff">
          <tr class="file-code-line">
            <td class="blob-line-nums">
              <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
<span id="L292" rel="#L292">292</span>
<span id="L293" rel="#L293">293</span>
<span id="L294" rel="#L294">294</span>
<span id="L295" rel="#L295">295</span>
<span id="L296" rel="#L296">296</span>
<span id="L297" rel="#L297">297</span>
<span id="L298" rel="#L298">298</span>
<span id="L299" rel="#L299">299</span>
<span id="L300" rel="#L300">300</span>
<span id="L301" rel="#L301">301</span>
<span id="L302" rel="#L302">302</span>
<span id="L303" rel="#L303">303</span>
<span id="L304" rel="#L304">304</span>
<span id="L305" rel="#L305">305</span>
<span id="L306" rel="#L306">306</span>
<span id="L307" rel="#L307">307</span>
<span id="L308" rel="#L308">308</span>
<span id="L309" rel="#L309">309</span>
<span id="L310" rel="#L310">310</span>
<span id="L311" rel="#L311">311</span>
<span id="L312" rel="#L312">312</span>
<span id="L313" rel="#L313">313</span>
<span id="L314" rel="#L314">314</span>
<span id="L315" rel="#L315">315</span>
<span id="L316" rel="#L316">316</span>
<span id="L317" rel="#L317">317</span>
<span id="L318" rel="#L318">318</span>
<span id="L319" rel="#L319">319</span>
<span id="L320" rel="#L320">320</span>
<span id="L321" rel="#L321">321</span>
<span id="L322" rel="#L322">322</span>
<span id="L323" rel="#L323">323</span>
<span id="L324" rel="#L324">324</span>
<span id="L325" rel="#L325">325</span>
<span id="L326" rel="#L326">326</span>
<span id="L327" rel="#L327">327</span>
<span id="L328" rel="#L328">328</span>
<span id="L329" rel="#L329">329</span>
<span id="L330" rel="#L330">330</span>
<span id="L331" rel="#L331">331</span>
<span id="L332" rel="#L332">332</span>
<span id="L333" rel="#L333">333</span>
<span id="L334" rel="#L334">334</span>
<span id="L335" rel="#L335">335</span>
<span id="L336" rel="#L336">336</span>
<span id="L337" rel="#L337">337</span>
<span id="L338" rel="#L338">338</span>
<span id="L339" rel="#L339">339</span>
<span id="L340" rel="#L340">340</span>
<span id="L341" rel="#L341">341</span>
<span id="L342" rel="#L342">342</span>
<span id="L343" rel="#L343">343</span>
<span id="L344" rel="#L344">344</span>
<span id="L345" rel="#L345">345</span>
<span id="L346" rel="#L346">346</span>
<span id="L347" rel="#L347">347</span>
<span id="L348" rel="#L348">348</span>
<span id="L349" rel="#L349">349</span>
<span id="L350" rel="#L350">350</span>
<span id="L351" rel="#L351">351</span>
<span id="L352" rel="#L352">352</span>
<span id="L353" rel="#L353">353</span>
<span id="L354" rel="#L354">354</span>
<span id="L355" rel="#L355">355</span>
<span id="L356" rel="#L356">356</span>
<span id="L357" rel="#L357">357</span>
<span id="L358" rel="#L358">358</span>
<span id="L359" rel="#L359">359</span>
<span id="L360" rel="#L360">360</span>
<span id="L361" rel="#L361">361</span>
<span id="L362" rel="#L362">362</span>
<span id="L363" rel="#L363">363</span>
<span id="L364" rel="#L364">364</span>
<span id="L365" rel="#L365">365</span>
<span id="L366" rel="#L366">366</span>
<span id="L367" rel="#L367">367</span>
<span id="L368" rel="#L368">368</span>
<span id="L369" rel="#L369">369</span>
<span id="L370" rel="#L370">370</span>
<span id="L371" rel="#L371">371</span>
<span id="L372" rel="#L372">372</span>
<span id="L373" rel="#L373">373</span>
<span id="L374" rel="#L374">374</span>
<span id="L375" rel="#L375">375</span>
<span id="L376" rel="#L376">376</span>
<span id="L377" rel="#L377">377</span>
<span id="L378" rel="#L378">378</span>
<span id="L379" rel="#L379">379</span>
<span id="L380" rel="#L380">380</span>
<span id="L381" rel="#L381">381</span>
<span id="L382" rel="#L382">382</span>
<span id="L383" rel="#L383">383</span>
<span id="L384" rel="#L384">384</span>
<span id="L385" rel="#L385">385</span>
<span id="L386" rel="#L386">386</span>
<span id="L387" rel="#L387">387</span>
<span id="L388" rel="#L388">388</span>
<span id="L389" rel="#L389">389</span>
<span id="L390" rel="#L390">390</span>
<span id="L391" rel="#L391">391</span>
<span id="L392" rel="#L392">392</span>
<span id="L393" rel="#L393">393</span>
<span id="L394" rel="#L394">394</span>
<span id="L395" rel="#L395">395</span>
<span id="L396" rel="#L396">396</span>
<span id="L397" rel="#L397">397</span>
<span id="L398" rel="#L398">398</span>
<span id="L399" rel="#L399">399</span>
<span id="L400" rel="#L400">400</span>
<span id="L401" rel="#L401">401</span>
<span id="L402" rel="#L402">402</span>
<span id="L403" rel="#L403">403</span>
<span id="L404" rel="#L404">404</span>
<span id="L405" rel="#L405">405</span>
<span id="L406" rel="#L406">406</span>
<span id="L407" rel="#L407">407</span>
<span id="L408" rel="#L408">408</span>
<span id="L409" rel="#L409">409</span>
<span id="L410" rel="#L410">410</span>
<span id="L411" rel="#L411">411</span>
<span id="L412" rel="#L412">412</span>
<span id="L413" rel="#L413">413</span>
<span id="L414" rel="#L414">414</span>
<span id="L415" rel="#L415">415</span>
<span id="L416" rel="#L416">416</span>
<span id="L417" rel="#L417">417</span>
<span id="L418" rel="#L418">418</span>
<span id="L419" rel="#L419">419</span>
<span id="L420" rel="#L420">420</span>
<span id="L421" rel="#L421">421</span>
<span id="L422" rel="#L422">422</span>
<span id="L423" rel="#L423">423</span>
<span id="L424" rel="#L424">424</span>
<span id="L425" rel="#L425">425</span>
<span id="L426" rel="#L426">426</span>
<span id="L427" rel="#L427">427</span>
<span id="L428" rel="#L428">428</span>
<span id="L429" rel="#L429">429</span>
<span id="L430" rel="#L430">430</span>
<span id="L431" rel="#L431">431</span>
<span id="L432" rel="#L432">432</span>
<span id="L433" rel="#L433">433</span>
<span id="L434" rel="#L434">434</span>
<span id="L435" rel="#L435">435</span>
<span id="L436" rel="#L436">436</span>
<span id="L437" rel="#L437">437</span>
<span id="L438" rel="#L438">438</span>
<span id="L439" rel="#L439">439</span>
<span id="L440" rel="#L440">440</span>
<span id="L441" rel="#L441">441</span>
<span id="L442" rel="#L442">442</span>
<span id="L443" rel="#L443">443</span>
<span id="L444" rel="#L444">444</span>
<span id="L445" rel="#L445">445</span>
<span id="L446" rel="#L446">446</span>
<span id="L447" rel="#L447">447</span>
<span id="L448" rel="#L448">448</span>
<span id="L449" rel="#L449">449</span>
<span id="L450" rel="#L450">450</span>
<span id="L451" rel="#L451">451</span>
<span id="L452" rel="#L452">452</span>
<span id="L453" rel="#L453">453</span>
<span id="L454" rel="#L454">454</span>
<span id="L455" rel="#L455">455</span>
<span id="L456" rel="#L456">456</span>
<span id="L457" rel="#L457">457</span>
<span id="L458" rel="#L458">458</span>
<span id="L459" rel="#L459">459</span>
<span id="L460" rel="#L460">460</span>
<span id="L461" rel="#L461">461</span>
<span id="L462" rel="#L462">462</span>
<span id="L463" rel="#L463">463</span>
<span id="L464" rel="#L464">464</span>
<span id="L465" rel="#L465">465</span>
<span id="L466" rel="#L466">466</span>
<span id="L467" rel="#L467">467</span>
<span id="L468" rel="#L468">468</span>
<span id="L469" rel="#L469">469</span>
<span id="L470" rel="#L470">470</span>
<span id="L471" rel="#L471">471</span>
<span id="L472" rel="#L472">472</span>
<span id="L473" rel="#L473">473</span>
<span id="L474" rel="#L474">474</span>
<span id="L475" rel="#L475">475</span>
<span id="L476" rel="#L476">476</span>
<span id="L477" rel="#L477">477</span>
<span id="L478" rel="#L478">478</span>
<span id="L479" rel="#L479">479</span>
<span id="L480" rel="#L480">480</span>
<span id="L481" rel="#L481">481</span>
<span id="L482" rel="#L482">482</span>
<span id="L483" rel="#L483">483</span>
<span id="L484" rel="#L484">484</span>
<span id="L485" rel="#L485">485</span>
<span id="L486" rel="#L486">486</span>
<span id="L487" rel="#L487">487</span>
<span id="L488" rel="#L488">488</span>
<span id="L489" rel="#L489">489</span>
<span id="L490" rel="#L490">490</span>
<span id="L491" rel="#L491">491</span>
<span id="L492" rel="#L492">492</span>
<span id="L493" rel="#L493">493</span>
<span id="L494" rel="#L494">494</span>
<span id="L495" rel="#L495">495</span>
<span id="L496" rel="#L496">496</span>
<span id="L497" rel="#L497">497</span>
<span id="L498" rel="#L498">498</span>
<span id="L499" rel="#L499">499</span>
<span id="L500" rel="#L500">500</span>
<span id="L501" rel="#L501">501</span>
<span id="L502" rel="#L502">502</span>
<span id="L503" rel="#L503">503</span>
<span id="L504" rel="#L504">504</span>
<span id="L505" rel="#L505">505</span>
<span id="L506" rel="#L506">506</span>
<span id="L507" rel="#L507">507</span>
<span id="L508" rel="#L508">508</span>
<span id="L509" rel="#L509">509</span>
<span id="L510" rel="#L510">510</span>
<span id="L511" rel="#L511">511</span>
<span id="L512" rel="#L512">512</span>
<span id="L513" rel="#L513">513</span>
<span id="L514" rel="#L514">514</span>
<span id="L515" rel="#L515">515</span>
<span id="L516" rel="#L516">516</span>
<span id="L517" rel="#L517">517</span>
<span id="L518" rel="#L518">518</span>
<span id="L519" rel="#L519">519</span>
<span id="L520" rel="#L520">520</span>
<span id="L521" rel="#L521">521</span>
<span id="L522" rel="#L522">522</span>
<span id="L523" rel="#L523">523</span>
<span id="L524" rel="#L524">524</span>
<span id="L525" rel="#L525">525</span>
<span id="L526" rel="#L526">526</span>
<span id="L527" rel="#L527">527</span>
<span id="L528" rel="#L528">528</span>
<span id="L529" rel="#L529">529</span>
<span id="L530" rel="#L530">530</span>
<span id="L531" rel="#L531">531</span>
<span id="L532" rel="#L532">532</span>

            </td>
            <td class="blob-line-code">
                    <div class="highlight"><pre><div class='line' id='LC1'>&nbsp;<span class="cm">/***************************************************************************** </span></div><div class='line' id='LC2'><span class="cm">  * Project: RooFit                                                           * </span></div><div class='line' id='LC3'><span class="cm">  *                                                                           * </span></div><div class='line' id='LC4'><span class="cm">  * This code was autogenerated by RooClassFactory                            * </span></div><div class='line' id='LC5'><span class="cm">  *****************************************************************************/</span> </div><div class='line' id='LC6'><br/></div><div class='line' id='LC7'>&nbsp;<span class="c1">// Your description goes here... </span></div><div class='line' id='LC8'><br/></div><div class='line' id='LC9'>&nbsp;<span class="err">#</span><span class="n">include</span> <span class="s">&quot;Riostream.h&quot;</span> </div><div class='line' id='LC10'><br/></div><div class='line' id='LC11'>&nbsp;<span class="err">#</span><span class="n">include</span> <span class="s">&quot;RooSpinZero_5D_ZH.h&quot;</span> </div><div class='line' id='LC12'>&nbsp;<span class="err">#</span><span class="n">include</span> <span class="s">&quot;RooAbsReal.h&quot;</span> </div><div class='line' id='LC13'>&nbsp;<span class="err">#</span><span class="n">include</span> <span class="s">&quot;RooAbsCategory.h&quot;</span> </div><div class='line' id='LC14'><br/></div><div class='line' id='LC15'><span class="cp">#include &lt;math.h&gt;</span></div><div class='line' id='LC16'><span class="cp">#include &quot;TMath.h&quot;</span></div><div class='line' id='LC17'><br/></div><div class='line' id='LC18'><span class="k">using</span> <span class="k">namespace</span> <span class="n">TMath</span><span class="p">;</span></div><div class='line' id='LC19'><br/></div><div class='line' id='LC20'><span class="k">enum</span> <span class="n">parameterizationList</span> <span class="p">{</span><span class="n">kMagPhase_As</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">kRealImag_Gs</span><span class="o">=</span><span class="mi">1</span><span class="p">,</span><span class="n">kFracPhase_Gs</span><span class="o">=</span><span class="mi">2</span><span class="p">,</span><span class="n">kNUMparameterizations</span><span class="o">=</span><span class="mi">3</span><span class="p">};</span></div><div class='line' id='LC21'><br/></div><div class='line' id='LC22'><br/></div><div class='line' id='LC23'>&nbsp;<span class="n">ClassImp</span><span class="p">(</span><span class="n">RooSpinZero_5D_ZH</span><span class="p">)</span> </div><div class='line' id='LC24'><br/></div><div class='line' id='LC25'>&nbsp;<span class="n">RooSpinZero_5D_ZH</span><span class="o">::</span><span class="n">RooSpinZero_5D_ZH</span><span class="p">(</span><span class="k">const</span> <span class="kt">char</span> <span class="o">*</span><span class="n">name</span><span class="p">,</span> <span class="k">const</span> <span class="kt">char</span> <span class="o">*</span><span class="n">title</span><span class="p">,</span> </div><div class='line' id='LC26'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_h1</span><span class="p">,</span></div><div class='line' id='LC27'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_h2</span><span class="p">,</span></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_hs</span><span class="p">,</span></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_Phi</span><span class="p">,</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_Phi1</span><span class="p">,</span></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_sqrts</span><span class="p">,</span></div><div class='line' id='LC32'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_mX</span><span class="p">,</span></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_mZ</span><span class="p">,</span></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_R1Val</span><span class="p">,</span></div><div class='line' id='LC35'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_R2Val</span><span class="p">,</span></div><div class='line' id='LC36'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="kt">int</span> <span class="n">_parameterization</span><span class="p">,</span></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_a1Val</span><span class="p">,</span></div><div class='line' id='LC38'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_phi1Val</span><span class="p">,</span></div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_a2Val</span><span class="p">,</span></div><div class='line' id='LC40'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_phi2Val</span><span class="p">,</span></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_a3Val</span><span class="p">,</span></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_phi3Val</span><span class="p">,</span></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g1Val</span><span class="p">,</span></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g2Val</span><span class="p">,</span></div><div class='line' id='LC45'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g3Val</span><span class="p">,</span></div><div class='line' id='LC46'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g4Val</span><span class="p">,</span></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g1ValIm</span><span class="p">,</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g2ValIm</span><span class="p">,</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g3ValIm</span><span class="p">,</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_g4ValIm</span><span class="p">,</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_fa2</span><span class="p">,</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_fa3</span><span class="p">,</span></div><div class='line' id='LC53'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_phia2</span><span class="p">,</span></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      <span class="n">RooAbsReal</span><span class="o">&amp;</span> <span class="n">_phia3</span> <span class="p">)</span><span class="o">:</span></div><div class='line' id='LC55'>&nbsp;&nbsp;&nbsp;<span class="n">RooAbsPdf</span><span class="p">(</span><span class="n">name</span><span class="p">,</span><span class="n">title</span><span class="p">),</span> </div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;<span class="n">h1</span><span class="p">(</span><span class="s">&quot;h1&quot;</span><span class="p">,</span><span class="s">&quot;h1&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_h1</span><span class="p">),</span></div><div class='line' id='LC57'>&nbsp;&nbsp;&nbsp;<span class="n">h2</span><span class="p">(</span><span class="s">&quot;h2&quot;</span><span class="p">,</span><span class="s">&quot;h2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_h2</span><span class="p">),</span></div><div class='line' id='LC58'>&nbsp;&nbsp;&nbsp;<span class="n">hs</span><span class="p">(</span><span class="s">&quot;hs&quot;</span><span class="p">,</span><span class="s">&quot;hs&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_hs</span><span class="p">),</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;<span class="n">Phi</span><span class="p">(</span><span class="s">&quot;Phi&quot;</span><span class="p">,</span><span class="s">&quot;Phi&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_Phi</span><span class="p">),</span></div><div class='line' id='LC60'>&nbsp;&nbsp;&nbsp;<span class="n">Phi1</span><span class="p">(</span><span class="s">&quot;Phi1&quot;</span><span class="p">,</span><span class="s">&quot;Phi1&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_Phi1</span><span class="p">),</span></div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;<span class="n">sqrts</span><span class="p">(</span><span class="s">&quot;sqrts&quot;</span><span class="p">,</span><span class="s">&quot;sqrts&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_sqrts</span><span class="p">),</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;<span class="n">mX</span><span class="p">(</span><span class="s">&quot;mX&quot;</span><span class="p">,</span><span class="s">&quot;mX&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_mX</span><span class="p">),</span></div><div class='line' id='LC63'>&nbsp;&nbsp;&nbsp;<span class="n">mZ</span><span class="p">(</span><span class="s">&quot;mZ&quot;</span><span class="p">,</span><span class="s">&quot;mZ&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_mZ</span><span class="p">),</span></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;<span class="n">R1Val</span><span class="p">(</span><span class="s">&quot;R1Val&quot;</span><span class="p">,</span><span class="s">&quot;R1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_R1Val</span><span class="p">),</span></div><div class='line' id='LC65'>&nbsp;&nbsp;&nbsp;<span class="n">R2Val</span><span class="p">(</span><span class="s">&quot;R2Val&quot;</span><span class="p">,</span><span class="s">&quot;R2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_R2Val</span><span class="p">),</span></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;<span class="n">parameterization</span><span class="p">(</span><span class="n">_parameterization</span><span class="p">),</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;<span class="n">a1Val</span><span class="p">(</span><span class="s">&quot;a1Val&quot;</span><span class="p">,</span><span class="s">&quot;a1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_a1Val</span><span class="p">),</span></div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;<span class="n">phi1Val</span><span class="p">(</span><span class="s">&quot;phi1Val&quot;</span><span class="p">,</span><span class="s">&quot;phi1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_phi1Val</span><span class="p">),</span></div><div class='line' id='LC69'>&nbsp;&nbsp;&nbsp;<span class="n">a2Val</span><span class="p">(</span><span class="s">&quot;a2Val&quot;</span><span class="p">,</span><span class="s">&quot;a2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_a2Val</span><span class="p">),</span></div><div class='line' id='LC70'>&nbsp;&nbsp;&nbsp;<span class="n">phi2Val</span><span class="p">(</span><span class="s">&quot;phi2Val&quot;</span><span class="p">,</span><span class="s">&quot;phi2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_phi2Val</span><span class="p">),</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;<span class="n">a3Val</span><span class="p">(</span><span class="s">&quot;a3Val&quot;</span><span class="p">,</span><span class="s">&quot;a3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_a3Val</span><span class="p">),</span></div><div class='line' id='LC72'>&nbsp;&nbsp;&nbsp;<span class="n">phi3Val</span><span class="p">(</span><span class="s">&quot;phi3Val&quot;</span><span class="p">,</span><span class="s">&quot;phi3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_phi3Val</span><span class="p">),</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;<span class="n">g1Val</span><span class="p">(</span><span class="s">&quot;g1Val&quot;</span><span class="p">,</span><span class="s">&quot;g1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g1Val</span><span class="p">),</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;<span class="n">g2Val</span><span class="p">(</span><span class="s">&quot;g2Val&quot;</span><span class="p">,</span><span class="s">&quot;g2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g2Val</span><span class="p">),</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;<span class="n">g3Val</span><span class="p">(</span><span class="s">&quot;g3Val&quot;</span><span class="p">,</span><span class="s">&quot;g3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g3Val</span><span class="p">),</span></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;<span class="n">g4Val</span><span class="p">(</span><span class="s">&quot;g4Val&quot;</span><span class="p">,</span><span class="s">&quot;g4Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g4Val</span><span class="p">),</span></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;<span class="n">g1ValIm</span><span class="p">(</span><span class="s">&quot;g1ValIm&quot;</span><span class="p">,</span><span class="s">&quot;g1ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g1ValIm</span><span class="p">),</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;<span class="n">g2ValIm</span><span class="p">(</span><span class="s">&quot;g2ValIm&quot;</span><span class="p">,</span><span class="s">&quot;g2ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g2ValIm</span><span class="p">),</span></div><div class='line' id='LC79'>&nbsp;&nbsp;&nbsp;<span class="n">g3ValIm</span><span class="p">(</span><span class="s">&quot;g3ValIm&quot;</span><span class="p">,</span><span class="s">&quot;g3ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g3ValIm</span><span class="p">),</span></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;<span class="n">g4ValIm</span><span class="p">(</span><span class="s">&quot;g4ValIm&quot;</span><span class="p">,</span><span class="s">&quot;g4ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_g4ValIm</span><span class="p">),</span></div><div class='line' id='LC81'>&nbsp;&nbsp;&nbsp;<span class="n">fa2</span><span class="p">(</span><span class="s">&quot;fa2&quot;</span><span class="p">,</span><span class="s">&quot;fa2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_fa2</span><span class="p">),</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;<span class="n">fa3</span><span class="p">(</span><span class="s">&quot;fa3&quot;</span><span class="p">,</span><span class="s">&quot;fa3&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_fa3</span><span class="p">),</span></div><div class='line' id='LC83'>&nbsp;&nbsp;&nbsp;<span class="n">phia2</span><span class="p">(</span><span class="s">&quot;phia2&quot;</span><span class="p">,</span><span class="s">&quot;phia2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_phia2</span><span class="p">),</span></div><div class='line' id='LC84'>&nbsp;&nbsp;&nbsp;<span class="n">phia3</span><span class="p">(</span><span class="s">&quot;phia3&quot;</span><span class="p">,</span><span class="s">&quot;phia3&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">_phia3</span><span class="p">)</span></div><div class='line' id='LC85'>&nbsp;<span class="p">{</span> </div><div class='line' id='LC86'>&nbsp;<span class="p">}</span> </div><div class='line' id='LC87'><br/></div><div class='line' id='LC88'><br/></div><div class='line' id='LC89'>&nbsp;<span class="n">RooSpinZero_5D_ZH</span><span class="o">::</span><span class="n">RooSpinZero_5D_ZH</span><span class="p">(</span><span class="k">const</span> <span class="n">RooSpinZero_5D_ZH</span><span class="o">&amp;</span> <span class="n">other</span><span class="p">,</span> <span class="k">const</span> <span class="kt">char</span><span class="o">*</span> <span class="n">name</span><span class="p">)</span> <span class="o">:</span>  </div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;<span class="n">RooAbsPdf</span><span class="p">(</span><span class="n">other</span><span class="p">,</span><span class="n">name</span><span class="p">),</span> </div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;<span class="n">h1</span><span class="p">(</span><span class="s">&quot;h1&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">h1</span><span class="p">),</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;<span class="n">h2</span><span class="p">(</span><span class="s">&quot;h2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">h2</span><span class="p">),</span></div><div class='line' id='LC93'>&nbsp;&nbsp;&nbsp;<span class="n">hs</span><span class="p">(</span><span class="s">&quot;hs&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">hs</span><span class="p">),</span></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;<span class="n">Phi</span><span class="p">(</span><span class="s">&quot;Phi&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">Phi</span><span class="p">),</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;<span class="n">Phi1</span><span class="p">(</span><span class="s">&quot;Phi1&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">Phi1</span><span class="p">),</span></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;<span class="n">sqrts</span><span class="p">(</span><span class="s">&quot;sqrts&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">sqrts</span><span class="p">),</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;<span class="n">mX</span><span class="p">(</span><span class="s">&quot;mX&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">mX</span><span class="p">),</span></div><div class='line' id='LC98'>&nbsp;&nbsp;&nbsp;<span class="n">mZ</span><span class="p">(</span><span class="s">&quot;mZ&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">mZ</span><span class="p">),</span></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;<span class="n">R1Val</span><span class="p">(</span><span class="s">&quot;R1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">R1Val</span><span class="p">),</span></div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;<span class="n">R2Val</span><span class="p">(</span><span class="s">&quot;R2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">R2Val</span><span class="p">),</span></div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;<span class="n">parameterization</span><span class="p">(</span><span class="n">other</span><span class="p">.</span><span class="n">parameterization</span><span class="p">),</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;<span class="n">a1Val</span><span class="p">(</span><span class="s">&quot;a1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">a1Val</span><span class="p">),</span></div><div class='line' id='LC103'>&nbsp;&nbsp;&nbsp;<span class="n">phi1Val</span><span class="p">(</span><span class="s">&quot;phi1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">phi1Val</span><span class="p">),</span></div><div class='line' id='LC104'>&nbsp;&nbsp;&nbsp;<span class="n">a2Val</span><span class="p">(</span><span class="s">&quot;a2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">a2Val</span><span class="p">),</span></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;<span class="n">phi2Val</span><span class="p">(</span><span class="s">&quot;phi2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">phi2Val</span><span class="p">),</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;<span class="n">a3Val</span><span class="p">(</span><span class="s">&quot;a3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">a3Val</span><span class="p">),</span></div><div class='line' id='LC107'>&nbsp;&nbsp;&nbsp;<span class="n">phi3Val</span><span class="p">(</span><span class="s">&quot;phi3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">phi3Val</span><span class="p">),</span></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;<span class="n">g1Val</span><span class="p">(</span><span class="s">&quot;g1Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g1Val</span><span class="p">),</span></div><div class='line' id='LC109'>&nbsp;&nbsp;&nbsp;<span class="n">g2Val</span><span class="p">(</span><span class="s">&quot;a2Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g2Val</span><span class="p">),</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;<span class="n">g3Val</span><span class="p">(</span><span class="s">&quot;g3Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g3Val</span><span class="p">),</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;<span class="n">g4Val</span><span class="p">(</span><span class="s">&quot;g4Val&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g4Val</span><span class="p">),</span></div><div class='line' id='LC112'>&nbsp;&nbsp;&nbsp;<span class="n">g1ValIm</span><span class="p">(</span><span class="s">&quot;g1ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g1ValIm</span><span class="p">),</span></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;<span class="n">g2ValIm</span><span class="p">(</span><span class="s">&quot;a2ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g2ValIm</span><span class="p">),</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;<span class="n">g3ValIm</span><span class="p">(</span><span class="s">&quot;g3ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g3ValIm</span><span class="p">),</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;<span class="n">g4ValIm</span><span class="p">(</span><span class="s">&quot;g4ValIm&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">g4ValIm</span><span class="p">),</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;<span class="n">fa2</span><span class="p">(</span><span class="s">&quot;fa2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">fa2</span><span class="p">),</span></div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;<span class="n">fa3</span><span class="p">(</span><span class="s">&quot;fa3&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">fa3</span><span class="p">),</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;<span class="n">phia2</span><span class="p">(</span><span class="s">&quot;phia2&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">phia2</span><span class="p">),</span></div><div class='line' id='LC119'>&nbsp;&nbsp;&nbsp;<span class="n">phia3</span><span class="p">(</span><span class="s">&quot;phia3&quot;</span><span class="p">,</span><span class="k">this</span><span class="p">,</span><span class="n">other</span><span class="p">.</span><span class="n">phia3</span><span class="p">)</span></div><div class='line' id='LC120'>&nbsp;<span class="p">{</span> </div><div class='line' id='LC121'>&nbsp;<span class="p">}</span> </div><div class='line' id='LC122'><br/></div><div class='line' id='LC123'><br/></div><div class='line' id='LC124'><br/></div><div class='line' id='LC125'>&nbsp;<span class="n">Double_t</span> <span class="n">RooSpinZero_5D_ZH</span><span class="o">::</span><span class="n">evaluate</span><span class="p">()</span> <span class="k">const</span> </div><div class='line' id='LC126'>&nbsp;<span class="p">{</span> </div><div class='line' id='LC127'><br/></div><div class='line' id='LC128'>&nbsp;&nbsp;&nbsp;<span class="c1">// these amplitudes are calculated based on comparing equations to the PRD paper</span></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;<span class="c1">// the beta and gamma are velocities of the Z in the C.O.M</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;<span class="c1">// Double_t betaValSquared =(pow(sqrts,2)-(pow(mX+mZ,2)))*(pow(sqrts,2)-(pow(mX-mZ,2)))/pow(pow(sqrts,2)-mX*mX+mZ*mZ,2);</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;<span class="c1">// Double_t gamma = 1./sqrt(1-betaValSquared);</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// http://prd.aps.org/pdf/PRD/v49/i1/p79_1</span></div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;<span class="c1">// Double_t f00 = gamma*gamma/(2+gamma*gamma);</span></div><div class='line' id='LC134'><br/></div><div class='line' id='LC135'><br/></div><div class='line' id='LC136'>&nbsp;&nbsp;&nbsp;<span class="c1">// below calcualtions are based on the H-&gt;ZZ amplitudes </span></div><div class='line' id='LC137'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">s</span><span class="o">=-</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="o">-</span><span class="n">sqrts</span><span class="o">*</span><span class="n">sqrts</span><span class="o">-</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="p">)</span><span class="o">/</span><span class="mf">2.</span><span class="p">;</span></div><div class='line' id='LC138'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">kappa</span><span class="o">=</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="mi">1000</span><span class="o">*</span><span class="mi">1000</span><span class="p">);</span></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">a1</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a2</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a3</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a1Im</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a2Im</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a3Im</span><span class="o">=</span><span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC141'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="nf">g1</span><span class="p">(</span><span class="mf">1.0</span><span class="p">),</span> <span class="n">g1Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g2</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g2Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g3</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g3Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g4</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g4Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">);</span></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC143'>&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kMagPhase_As</span><span class="p">){</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1</span><span class="o">=</span><span class="n">a1Val</span><span class="p">;</span></div><div class='line' id='LC145'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1Im</span><span class="o">=</span><span class="n">phi1Val</span><span class="p">;</span></div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2</span><span class="o">=</span><span class="n">a2Val</span><span class="p">;</span></div><div class='line' id='LC147'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2Im</span><span class="o">=</span><span class="n">phi2Val</span><span class="p">;</span></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3</span><span class="o">=</span><span class="n">a3Val</span><span class="p">;</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3Im</span><span class="o">=</span><span class="n">phi3Val</span><span class="p">;</span></div><div class='line' id='LC150'>&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC151'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kFracPhase_Gs</span><span class="p">){</span></div><div class='line' id='LC153'><br/></div><div class='line' id='LC154'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">nanval</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">fa2</span> <span class="o">-</span> <span class="n">fa3</span><span class="p">);</span></div><div class='line' id='LC155'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">nanval</span> <span class="o">!=</span> <span class="n">nanval</span><span class="p">)</span> <span class="k">return</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC156'><br/></div><div class='line' id='LC157'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// convert fraction and phase to g1,g2...etc</span></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC numbers at 250 GeV at mH= 125 GeV (narrow Z width approximation)</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">0.981396</span><span class="p">;</span> <span class="c1">// was 0.94696 at 126 GeV</span></div><div class='line' id='LC160'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">33.4674</span><span class="p">;</span>  <span class="c1">// was 32.1981 at 126 GeV</span></div><div class='line' id='LC161'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">7.9229</span><span class="p">;</span>   <span class="c1">// was 7.45502 at 126 GeV</span></div><div class='line' id='LC162'><br/></div><div class='line' id='LC163'><br/></div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 350 GeV at mH = 125 GeV</span></div><div class='line' id='LC165'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">350.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC166'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">1.48872</span><span class="p">;</span> </div><div class='line' id='LC167'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">125.387</span><span class="p">;</span>  </div><div class='line' id='LC168'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">75.3199</span><span class="p">;</span></div><div class='line' id='LC169'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC170'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC171'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 500 GeV at mH = 125 GeV</span></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">500.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC173'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">2.57246</span><span class="p">;</span> </div><div class='line' id='LC174'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">516.556</span><span class="p">;</span>  </div><div class='line' id='LC175'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">414.378</span><span class="p">;</span></div><div class='line' id='LC176'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC178'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 1000 GeV at mH = 125 GeV</span></div><div class='line' id='LC179'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">1000.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">8.95721</span><span class="p">;</span> </div><div class='line' id='LC181'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">8208.91</span><span class="p">;</span>  </div><div class='line' id='LC182'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">7800.2</span><span class="p">;</span></div><div class='line' id='LC183'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC184'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC186'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC188'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC190'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">fa2</span> <span class="o">==</span> <span class="mf">1.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC191'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC192'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span></div><div class='line' id='LC193'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC194'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="k">if</span> <span class="p">(</span> <span class="n">fa3</span> <span class="o">==</span> <span class="mf">1.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC196'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span> </div><div class='line' id='LC198'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC199'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Mag</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="n">fa2</span><span class="o">/</span><span class="p">(</span><span class="mf">1.</span><span class="o">-</span><span class="n">fa2</span><span class="o">-</span><span class="n">fa3</span><span class="p">))</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">sigma1_e</span><span class="o">/</span><span class="n">sigma2_e</span><span class="p">);</span> </div><div class='line' id='LC200'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Mag</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="n">fa3</span><span class="o">/</span><span class="p">(</span><span class="mf">1.</span><span class="o">-</span><span class="n">fa2</span><span class="o">-</span><span class="n">fa3</span><span class="p">))</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">sigma1_e</span><span class="o">/</span><span class="n">sigma4_e</span><span class="p">);</span> </div><div class='line' id='LC201'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1</span>   <span class="o">=</span> <span class="n">g1Mag</span><span class="p">;</span></div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Im</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC205'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g2Mag</span><span class="o">*</span><span class="n">cos</span><span class="p">(</span><span class="n">phia2</span><span class="p">);</span></div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g2Mag</span><span class="o">*</span><span class="n">sin</span><span class="p">(</span><span class="n">phia2</span><span class="p">);</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3</span>   <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC208'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3Im</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC209'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g4Mag</span><span class="o">*</span><span class="n">cos</span><span class="p">(</span><span class="n">phia3</span><span class="p">);</span></div><div class='line' id='LC210'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g4Mag</span><span class="o">*</span><span class="n">sin</span><span class="p">(</span><span class="n">phia3</span><span class="p">);</span></div><div class='line' id='LC211'><br/></div><div class='line' id='LC212'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span><span class="k">else</span> <span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kRealImag_Gs</span><span class="p">){</span></div><div class='line' id='LC213'><br/></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1</span>   <span class="o">=</span>  <span class="n">g1Val</span><span class="p">;</span></div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Im</span> <span class="o">=</span>  <span class="n">g1ValIm</span><span class="p">;</span></div><div class='line' id='LC216'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2</span>   <span class="o">=</span>  <span class="o">-</span> <span class="n">g2Val</span><span class="p">;</span></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Im</span> <span class="o">=</span>  <span class="o">-</span> <span class="n">g2ValIm</span><span class="p">;</span></div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3</span>   <span class="o">=</span> <span class="n">g3Val</span><span class="p">;</span></div><div class='line' id='LC219'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3Im</span> <span class="o">=</span>  <span class="n">g3ValIm</span><span class="p">;</span></div><div class='line' id='LC220'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g4Val</span><span class="p">;</span></div><div class='line' id='LC221'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g4ValIm</span><span class="p">;</span></div><div class='line' id='LC222'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC224'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1</span> <span class="o">=</span> <span class="n">g1</span><span class="o">*</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g2</span><span class="o">*</span><span class="mf">2.</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g3</span><span class="o">*</span><span class="n">kappa</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">);</span></div><div class='line' id='LC225'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1Im</span> <span class="o">=</span> <span class="n">g1Im</span><span class="o">*</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g2Im</span><span class="o">*</span><span class="mf">2.</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g3Im</span><span class="o">*</span><span class="n">kappa</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">);</span></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g2</span> <span class="o">-</span> <span class="n">g3</span><span class="o">*</span><span class="n">kappa</span><span class="p">;</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2Im</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g2Im</span> <span class="o">-</span> <span class="n">g3Im</span><span class="o">*</span><span class="n">kappa</span><span class="p">;</span></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g4</span><span class="p">;</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3Im</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g4Im</span><span class="p">;</span></div><div class='line' id='LC230'><br/></div><div class='line' id='LC231'>&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC232'>&nbsp;&nbsp;</div><div class='line' id='LC233'><br/></div><div class='line' id='LC234'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">x</span> <span class="o">=</span> <span class="n">pow</span><span class="p">((</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="o">-</span><span class="n">sqrts</span><span class="o">*</span><span class="n">sqrts</span><span class="o">-</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="mf">2.</span><span class="o">*</span><span class="n">sqrts</span><span class="o">*</span><span class="n">mZ</span><span class="p">),</span><span class="mi">2</span><span class="p">)</span><span class="o">-</span><span class="mi">1</span><span class="p">;</span></div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">A00Real</span> <span class="o">=</span> <span class="o">-</span> <span class="p">(</span><span class="n">a1</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">x</span><span class="p">)</span> <span class="o">+</span> <span class="n">a2</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">A00Imag</span> <span class="o">=</span> <span class="o">-</span> <span class="p">(</span><span class="n">a1Im</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">x</span><span class="p">)</span> <span class="o">+</span> <span class="n">a2Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC238'><br/></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Ap0Real</span> <span class="o">=</span> <span class="n">a1</span> <span class="o">-</span> <span class="n">a3Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Ap0Imag</span> <span class="o">=</span> <span class="n">a1Im</span> <span class="o">+</span> <span class="n">a3</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC241'><br/></div><div class='line' id='LC242'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Am0Real</span> <span class="o">=</span> <span class="n">a1</span> <span class="o">+</span> <span class="n">a3Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Am0Imag</span> <span class="o">=</span> <span class="n">a1Im</span> <span class="o">-</span> <span class="n">a3</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC244'><br/></div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">f00</span> <span class="o">=</span> <span class="n">A00Real</span><span class="o">*</span><span class="n">A00Real</span> <span class="o">+</span> <span class="n">A00Imag</span><span class="o">*</span><span class="n">A00Imag</span><span class="p">;</span></div><div class='line' id='LC246'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">fp0</span> <span class="o">=</span> <span class="n">Ap0Real</span><span class="o">*</span><span class="n">Ap0Real</span> <span class="o">+</span> <span class="n">Ap0Imag</span><span class="o">*</span><span class="n">Ap0Imag</span><span class="p">;</span></div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">fm0</span> <span class="o">=</span> <span class="n">Am0Real</span><span class="o">*</span><span class="n">Am0Real</span> <span class="o">+</span> <span class="n">Am0Imag</span><span class="o">*</span><span class="n">Am0Imag</span><span class="p">;</span></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">ftotal</span> <span class="o">=</span> <span class="n">f00</span> <span class="o">+</span> <span class="n">fp0</span> <span class="o">+</span> <span class="n">fm0</span><span class="p">;</span></div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;<span class="c1">// normalize to the total</span></div><div class='line' id='LC252'>&nbsp;&nbsp;&nbsp;<span class="n">f00</span> <span class="o">=</span> <span class="n">f00</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC253'>&nbsp;&nbsp;&nbsp;<span class="n">fp0</span> <span class="o">=</span> <span class="n">fp0</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;<span class="n">fm0</span> <span class="o">=</span> <span class="n">fm0</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC255'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC256'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phi00</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">A00Imag</span><span class="p">,</span><span class="n">A00Real</span><span class="p">);</span></div><div class='line' id='LC257'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phip0</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">Ap0Imag</span><span class="p">,</span><span class="n">Ap0Real</span><span class="p">)</span><span class="o">-</span><span class="n">phi00</span><span class="p">;</span></div><div class='line' id='LC258'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phim0</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">Am0Imag</span><span class="p">,</span><span class="n">Am0Real</span><span class="p">)</span><span class="o">-</span><span class="n">phi00</span><span class="p">;</span></div><div class='line' id='LC259'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC260'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">value</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC261'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC262'><br/></div><div class='line' id='LC263'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">f00</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">)))</span><span class="o">/</span><span class="mf">4.</span><span class="p">;</span></div><div class='line' id='LC264'><br/></div><div class='line' id='LC265'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">fp0</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">-</span> <span class="mi">2</span><span class="o">*</span><span class="n">h1</span><span class="o">*</span><span class="n">R1Val</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">2</span><span class="o">*</span><span class="n">h2</span><span class="o">*</span><span class="n">R2Val</span><span class="p">))</span><span class="o">/</span><span class="mf">16.</span><span class="p">;</span></div><div class='line' id='LC266'><br/></div><div class='line' id='LC267'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">fm0</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">2</span><span class="o">*</span><span class="n">h1</span><span class="o">*</span><span class="n">R1Val</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">-</span> <span class="mi">2</span><span class="o">*</span><span class="n">h2</span><span class="o">*</span><span class="n">R2Val</span><span class="p">))</span><span class="o">/</span><span class="mf">16.</span><span class="p">;</span></div><div class='line' id='LC268'><br/></div><div class='line' id='LC269'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="o">-</span><span class="p">(</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">f00</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fp0</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="p">(</span><span class="n">h1</span> <span class="o">-</span> <span class="n">R1Val</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="n">h2</span> <span class="o">+</span> <span class="n">R2Val</span><span class="p">)</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="n">Phi</span> <span class="o">+</span> <span class="n">phip0</span><span class="p">))</span><span class="o">/</span><span class="mf">4.</span><span class="p">;</span></div><div class='line' id='LC270'><br/></div><div class='line' id='LC271'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="o">-</span><span class="p">(</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">f00</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fm0</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="p">(</span><span class="n">h1</span> <span class="o">+</span> <span class="n">R1Val</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="n">h2</span> <span class="o">-</span> <span class="n">R2Val</span><span class="p">)</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="n">Phi</span> <span class="o">-</span> <span class="n">phim0</span><span class="p">))</span><span class="o">/</span><span class="mf">4.</span><span class="p">;</span></div><div class='line' id='LC272'><br/></div><div class='line' id='LC273'>&nbsp;&nbsp;&nbsp;<span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fm0</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fp0</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="mi">2</span><span class="o">*</span><span class="n">Phi</span> <span class="o">-</span> <span class="n">phim0</span> <span class="o">+</span> <span class="n">phip0</span><span class="p">))</span><span class="o">/</span><span class="mf">8.</span><span class="p">;</span></div><div class='line' id='LC274'><br/></div><div class='line' id='LC275'><br/></div><div class='line' id='LC276'>&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">value</span> <span class="p">;</span> </div><div class='line' id='LC277'>&nbsp;<span class="p">}</span> </div><div class='line' id='LC278'><br/></div><div class='line' id='LC279'><span class="n">Int_t</span> <span class="n">RooSpinZero_5D_ZH</span><span class="o">::</span><span class="n">getAnalyticalIntegral</span><span class="p">(</span><span class="n">RooArgSet</span><span class="o">&amp;</span> <span class="n">allVars</span><span class="p">,</span> <span class="n">RooArgSet</span><span class="o">&amp;</span> <span class="n">analVars</span><span class="p">,</span> <span class="k">const</span> <span class="kt">char</span><span class="o">*</span> <span class="cm">/*rangeName*/</span><span class="p">)</span> <span class="k">const</span></div><div class='line' id='LC280'><span class="p">{</span></div><div class='line' id='LC281'><br/></div><div class='line' id='LC282'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">RooArgSet</span><span class="p">(</span><span class="o">*</span><span class="n">h1</span><span class="p">.</span><span class="n">absArg</span><span class="p">(),</span><span class="o">*</span><span class="n">h2</span><span class="p">.</span><span class="n">absArg</span><span class="p">(),</span><span class="o">*</span><span class="n">hs</span><span class="p">.</span><span class="n">absArg</span><span class="p">(),</span><span class="o">*</span><span class="n">Phi1</span><span class="p">.</span><span class="n">absArg</span><span class="p">(),</span><span class="o">*</span><span class="n">Phi</span><span class="p">.</span><span class="n">absArg</span><span class="p">())))</span> <span class="k">return</span> <span class="mi">6</span> <span class="p">;</span></div><div class='line' id='LC283'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">h1</span><span class="p">,</span><span class="n">h2</span><span class="p">,</span><span class="n">Phi1</span><span class="p">,</span><span class="n">Phi</span><span class="p">))</span> <span class="k">return</span> <span class="mi">1</span> <span class="p">;</span></div><div class='line' id='LC284'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">h1</span><span class="p">,</span><span class="n">h2</span><span class="p">,</span><span class="n">hs</span><span class="p">,</span><span class="n">Phi1</span><span class="p">))</span> <span class="k">return</span> <span class="mi">5</span> <span class="p">;</span></div><div class='line' id='LC285'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">h1</span><span class="p">,</span><span class="n">hs</span><span class="p">,</span><span class="n">Phi1</span><span class="p">,</span><span class="n">Phi</span><span class="p">))</span> <span class="k">return</span> <span class="mi">3</span> <span class="p">;</span></div><div class='line' id='LC286'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">h2</span><span class="p">,</span><span class="n">hs</span><span class="p">,</span><span class="n">Phi1</span><span class="p">,</span><span class="n">Phi</span><span class="p">))</span> <span class="k">return</span> <span class="mi">4</span> <span class="p">;</span></div><div class='line' id='LC287'>&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">matchArgs</span><span class="p">(</span><span class="n">allVars</span><span class="p">,</span><span class="n">analVars</span><span class="p">,</span><span class="n">h1</span><span class="p">,</span><span class="n">h2</span><span class="p">,</span><span class="n">hs</span><span class="p">,</span><span class="n">Phi</span><span class="p">))</span> <span class="k">return</span> <span class="mi">2</span> <span class="p">;</span></div><div class='line' id='LC288'><br/></div><div class='line' id='LC289'>&nbsp;&nbsp;<span class="k">return</span> <span class="mi">0</span> <span class="p">;</span></div><div class='line' id='LC290'><span class="p">}</span></div><div class='line' id='LC291'><br/></div><div class='line' id='LC292'><span class="n">Double_t</span> <span class="n">RooSpinZero_5D_ZH</span><span class="o">::</span><span class="n">analyticalIntegral</span><span class="p">(</span><span class="n">Int_t</span> <span class="n">code</span><span class="p">,</span> <span class="k">const</span> <span class="kt">char</span><span class="o">*</span> <span class="n">rangeName</span><span class="p">)</span> <span class="k">const</span></div><div class='line' id='LC293'><span class="p">{</span></div><div class='line' id='LC294'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC295'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// the beta and gamma are velocities of the Z in the C.O.M</span></div><div class='line' id='LC296'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">betaValSquared</span> <span class="o">=</span><span class="p">(</span><span class="n">pow</span><span class="p">(</span><span class="n">sqrts</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span><span class="o">-</span><span class="p">(</span><span class="n">pow</span><span class="p">(</span><span class="n">mX</span><span class="o">+</span><span class="n">mZ</span><span class="p">,</span><span class="mi">2</span><span class="p">)))</span><span class="o">*</span><span class="p">(</span><span class="n">pow</span><span class="p">(</span><span class="n">sqrts</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span><span class="o">-</span><span class="p">(</span><span class="n">pow</span><span class="p">(</span><span class="n">mX</span><span class="o">-</span><span class="n">mZ</span><span class="p">,</span><span class="mi">2</span><span class="p">)))</span><span class="o">/</span><span class="n">pow</span><span class="p">(</span><span class="n">pow</span><span class="p">(</span><span class="n">sqrts</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span><span class="o">-</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="o">+</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="p">,</span><span class="mi">2</span><span class="p">);</span></div><div class='line' id='LC297'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">gamma</span> <span class="o">=</span> <span class="mf">1.</span><span class="o">/</span><span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span><span class="o">-</span><span class="n">betaValSquared</span><span class="p">);</span></div><div class='line' id='LC298'><br/></div><div class='line' id='LC299'>&nbsp;&nbsp;&nbsp;<span class="c1">// these amplitudes are calculated based on comparing equations to the PRD paper</span></div><div class='line' id='LC300'>&nbsp;&nbsp;&nbsp;<span class="c1">// http://prd.aps.org/pdf/PRD/v49/i1/p79_1</span></div><div class='line' id='LC301'>&nbsp;&nbsp;&nbsp;<span class="c1">// Double_t f00 = gamma*gamma/(2+gamma*gamma);</span></div><div class='line' id='LC302'>&nbsp;&nbsp;</div><div class='line' id='LC303'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// below calcualtions are based on the H-&gt;ZZ amplitudes </span></div><div class='line' id='LC304'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">s</span><span class="o">=-</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="o">-</span><span class="n">sqrts</span><span class="o">*</span><span class="n">sqrts</span><span class="o">-</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="p">)</span><span class="o">/</span><span class="mf">2.</span><span class="p">;</span></div><div class='line' id='LC305'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">kappa</span><span class="o">=</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="mi">1000</span><span class="o">*</span><span class="mi">1000</span><span class="p">);</span></div><div class='line' id='LC306'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC307'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">a1</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a2</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a3</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a1Im</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a2Im</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span><span class="n">a3Im</span><span class="o">=</span><span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC308'>&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="nf">g1</span><span class="p">(</span><span class="mf">1.0</span><span class="p">),</span> <span class="n">g1Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g2</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g2Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g3</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g3Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g4</span><span class="p">(</span><span class="mf">0.</span><span class="p">),</span> <span class="n">g4Im</span><span class="p">(</span><span class="mf">0.</span><span class="p">);</span></div><div class='line' id='LC309'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC310'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC311'>&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kMagPhase_As</span><span class="p">){</span></div><div class='line' id='LC312'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1</span><span class="o">=</span><span class="n">a1Val</span><span class="p">;</span></div><div class='line' id='LC313'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1Im</span><span class="o">=</span><span class="n">phi1Val</span><span class="p">;</span></div><div class='line' id='LC314'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2</span><span class="o">=</span><span class="n">a2Val</span><span class="p">;</span></div><div class='line' id='LC315'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2Im</span><span class="o">=</span><span class="n">phi2Val</span><span class="p">;</span></div><div class='line' id='LC316'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3</span><span class="o">=</span><span class="n">a3Val</span><span class="p">;</span></div><div class='line' id='LC317'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3Im</span><span class="o">=</span><span class="n">phi3Val</span><span class="p">;</span></div><div class='line' id='LC318'>&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC319'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC320'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kFracPhase_Gs</span><span class="p">){</span></div><div class='line' id='LC321'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// </span></div><div class='line' id='LC322'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// place holder, not correct</span></div><div class='line' id='LC323'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// </span></div><div class='line' id='LC324'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">nanval</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span> <span class="o">-</span> <span class="n">fa2</span> <span class="o">-</span> <span class="n">fa3</span><span class="p">);</span></div><div class='line' id='LC325'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">nanval</span> <span class="o">!=</span> <span class="n">nanval</span><span class="p">)</span> <span class="k">return</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC326'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC327'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// convert fraction and phase to g1,g2...etc</span></div><div class='line' id='LC328'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC numbers at 250 GeV at mH= 125 GeV (narrow Z width approximation)</span></div><div class='line' id='LC329'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">0.981396</span><span class="p">;</span> <span class="c1">// was 0.94696 at 126 GeV</span></div><div class='line' id='LC330'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">33.4674</span><span class="p">;</span>  <span class="c1">// was 32.1981 at 126 GeV</span></div><div class='line' id='LC331'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">7.9229</span><span class="p">;</span>   <span class="c1">// was 7.45502 at 126 GeV</span></div><div class='line' id='LC332'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC333'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 350 GeV at mH = 125 GeV</span></div><div class='line' id='LC334'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">350.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC335'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">1.48872</span><span class="p">;</span> </div><div class='line' id='LC336'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">125.387</span><span class="p">;</span>  </div><div class='line' id='LC337'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">75.3199</span><span class="p">;</span></div><div class='line' id='LC338'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC339'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC340'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 500 GeV at mH = 125 GeV</span></div><div class='line' id='LC341'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">500.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC342'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">2.57246</span><span class="p">;</span> </div><div class='line' id='LC343'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">516.556</span><span class="p">;</span>  </div><div class='line' id='LC344'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">414.378</span><span class="p">;</span></div><div class='line' id='LC345'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC346'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC347'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// ILC nubmers at 1000 GeV at mH = 125 GeV</span></div><div class='line' id='LC348'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">sqrts</span> <span class="o">==</span> <span class="mf">1000.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC349'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma1_e</span> <span class="o">=</span> <span class="mf">8.95721</span><span class="p">;</span> </div><div class='line' id='LC350'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma2_e</span> <span class="o">=</span> <span class="mf">8208.91</span><span class="p">;</span>  </div><div class='line' id='LC351'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">sigma4_e</span> <span class="o">=</span> <span class="mf">7800.2</span><span class="p">;</span></div><div class='line' id='LC352'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC353'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC354'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span></div><div class='line' id='LC355'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC356'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC357'><br/></div><div class='line' id='LC358'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">fa2</span> <span class="o">==</span> <span class="mf">1.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC359'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC360'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span></div><div class='line' id='LC361'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span> </div><div class='line' id='LC362'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="k">if</span> <span class="p">(</span> <span class="n">fa3</span> <span class="o">==</span> <span class="mf">1.</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC363'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g1Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC364'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g2Mag</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC365'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g4Mag</span> <span class="o">=</span> <span class="mf">1.</span><span class="p">;</span> </div><div class='line' id='LC366'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span> <span class="k">else</span> <span class="p">{</span></div><div class='line' id='LC367'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g2Mag</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="n">fa2</span><span class="o">/</span><span class="p">(</span><span class="mf">1.</span><span class="o">-</span><span class="n">fa2</span><span class="o">-</span><span class="n">fa3</span><span class="p">))</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">sigma1_e</span><span class="o">/</span><span class="n">sigma2_e</span><span class="p">);</span> </div><div class='line' id='LC368'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">g4Mag</span> <span class="o">=</span> <span class="n">sqrt</span><span class="p">(</span><span class="n">fa3</span><span class="o">/</span><span class="p">(</span><span class="mf">1.</span><span class="o">-</span><span class="n">fa2</span><span class="o">-</span><span class="n">fa3</span><span class="p">))</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">sigma1_e</span><span class="o">/</span><span class="n">sigma4_e</span><span class="p">);</span> </div><div class='line' id='LC369'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC370'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC371'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1</span>   <span class="o">=</span> <span class="n">g1Mag</span><span class="p">;</span></div><div class='line' id='LC372'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Im</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC373'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g2Mag</span><span class="o">*</span><span class="n">cos</span><span class="p">(</span><span class="n">phia2</span><span class="p">);</span></div><div class='line' id='LC374'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g2Mag</span><span class="o">*</span><span class="n">sin</span><span class="p">(</span><span class="n">phia2</span><span class="p">);</span></div><div class='line' id='LC375'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3</span>   <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC376'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3Im</span> <span class="o">=</span> <span class="mf">0.0</span><span class="p">;</span></div><div class='line' id='LC377'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g4Mag</span><span class="o">*</span><span class="n">cos</span><span class="p">(</span><span class="n">phia3</span><span class="p">);</span></div><div class='line' id='LC378'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g4Mag</span><span class="o">*</span><span class="n">sin</span><span class="p">(</span><span class="n">phia3</span><span class="p">);</span></div><div class='line' id='LC379'><br/></div><div class='line' id='LC380'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span><span class="k">else</span> <span class="k">if</span><span class="p">(</span><span class="n">parameterization</span><span class="o">==</span><span class="n">kRealImag_Gs</span><span class="p">){</span></div><div class='line' id='LC381'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC382'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1</span>   <span class="o">=</span>  <span class="n">g1Val</span><span class="p">;</span></div><div class='line' id='LC383'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g1Im</span> <span class="o">=</span>  <span class="n">g1ValIm</span><span class="p">;</span></div><div class='line' id='LC384'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2</span>   <span class="o">=</span>  <span class="o">-</span><span class="n">g2Val</span><span class="p">;</span></div><div class='line' id='LC385'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g2Im</span> <span class="o">=</span>  <span class="o">-</span><span class="n">g2ValIm</span><span class="p">;</span></div><div class='line' id='LC386'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3</span>   <span class="o">=</span> <span class="n">g3Val</span><span class="p">;</span></div><div class='line' id='LC387'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g3Im</span> <span class="o">=</span> <span class="n">g3ValIm</span><span class="p">;</span></div><div class='line' id='LC388'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4</span>   <span class="o">=</span> <span class="o">-</span> <span class="n">g4Val</span><span class="p">;</span></div><div class='line' id='LC389'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g4Im</span> <span class="o">=</span> <span class="o">-</span> <span class="n">g4ValIm</span><span class="p">;</span></div><div class='line' id='LC390'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC391'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC392'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1</span> <span class="o">=</span> <span class="n">g1</span><span class="o">*</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g2</span><span class="o">*</span><span class="mf">2.</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g3</span><span class="o">*</span><span class="n">kappa</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">);</span></div><div class='line' id='LC393'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a1Im</span> <span class="o">=</span> <span class="n">g1Im</span><span class="o">*</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g2Im</span><span class="o">*</span><span class="mf">2.</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span> <span class="o">+</span> <span class="n">g3Im</span><span class="o">*</span><span class="n">kappa</span><span class="o">*</span><span class="n">s</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">);</span></div><div class='line' id='LC394'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g2</span> <span class="o">-</span> <span class="n">g3</span><span class="o">*</span><span class="n">kappa</span><span class="p">;</span></div><div class='line' id='LC395'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a2Im</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g2Im</span> <span class="o">-</span> <span class="n">g3Im</span><span class="o">*</span><span class="n">kappa</span><span class="p">;</span></div><div class='line' id='LC396'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3</span> <span class="o">=</span>  <span class="o">-</span> <span class="mf">2.</span><span class="o">*</span><span class="n">g4</span><span class="p">;</span></div><div class='line' id='LC397'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">a3Im</span> <span class="o">=</span> <span class="o">-</span><span class="mf">2.</span><span class="o">*</span><span class="n">g4Im</span><span class="p">;</span></div><div class='line' id='LC398'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC399'>&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC400'>&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC401'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">x</span> <span class="o">=</span> <span class="n">pow</span><span class="p">((</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="o">-</span><span class="n">sqrts</span><span class="o">*</span><span class="n">sqrts</span><span class="o">-</span><span class="n">mZ</span><span class="o">*</span><span class="n">mZ</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="mf">2.</span><span class="o">*</span><span class="n">sqrts</span><span class="o">*</span><span class="n">mZ</span><span class="p">),</span><span class="mi">2</span><span class="p">)</span><span class="o">-</span><span class="mi">1</span><span class="p">;</span></div><div class='line' id='LC402'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC403'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">A00Real</span> <span class="o">=</span> <span class="o">-</span> <span class="p">(</span><span class="n">a1</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">x</span><span class="p">)</span> <span class="o">+</span> <span class="n">a2</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC404'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">A00Imag</span> <span class="o">=</span> <span class="o">-</span> <span class="p">(</span><span class="n">a1Im</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="mi">1</span><span class="o">+</span><span class="n">x</span><span class="p">)</span> <span class="o">+</span> <span class="n">a2Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC405'><br/></div><div class='line' id='LC406'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Ap0Real</span> <span class="o">=</span> <span class="n">a1</span> <span class="o">-</span> <span class="n">a3Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC407'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Ap0Imag</span> <span class="o">=</span> <span class="n">a1Im</span> <span class="o">+</span> <span class="n">a3</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC408'><br/></div><div class='line' id='LC409'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Am0Real</span> <span class="o">=</span> <span class="n">a1</span> <span class="o">+</span> <span class="n">a3Im</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC410'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">Am0Imag</span> <span class="o">=</span> <span class="n">a1Im</span> <span class="o">-</span> <span class="n">a3</span><span class="o">*</span><span class="p">(</span><span class="n">mZ</span><span class="o">*</span><span class="n">sqrts</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">mX</span><span class="o">*</span><span class="n">mX</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">x</span><span class="p">);</span></div><div class='line' id='LC411'><br/></div><div class='line' id='LC412'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">f00</span> <span class="o">=</span> <span class="n">A00Real</span><span class="o">*</span><span class="n">A00Real</span> <span class="o">+</span> <span class="n">A00Imag</span><span class="o">*</span><span class="n">A00Imag</span><span class="p">;</span></div><div class='line' id='LC413'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">fp0</span> <span class="o">=</span> <span class="n">Ap0Real</span><span class="o">*</span><span class="n">Ap0Real</span> <span class="o">+</span> <span class="n">Ap0Imag</span><span class="o">*</span><span class="n">Ap0Imag</span><span class="p">;</span></div><div class='line' id='LC414'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">fm0</span> <span class="o">=</span> <span class="n">Am0Real</span><span class="o">*</span><span class="n">Am0Real</span> <span class="o">+</span> <span class="n">Am0Imag</span><span class="o">*</span><span class="n">Am0Imag</span><span class="p">;</span></div><div class='line' id='LC415'><br/></div><div class='line' id='LC416'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">ftotal</span> <span class="o">=</span> <span class="n">f00</span> <span class="o">+</span> <span class="n">fp0</span> <span class="o">+</span> <span class="n">fm0</span><span class="p">;</span></div><div class='line' id='LC417'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC418'>&nbsp;&nbsp;&nbsp;<span class="c1">// normalize to the total</span></div><div class='line' id='LC419'>&nbsp;&nbsp;&nbsp;<span class="n">f00</span> <span class="o">=</span> <span class="n">f00</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC420'>&nbsp;&nbsp;&nbsp;<span class="n">fp0</span> <span class="o">=</span> <span class="n">fp0</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC421'>&nbsp;&nbsp;&nbsp;<span class="n">fm0</span> <span class="o">=</span> <span class="n">fm0</span> <span class="o">/</span> <span class="n">ftotal</span><span class="p">;</span></div><div class='line' id='LC422'>&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC423'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phi00</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">A00Imag</span><span class="p">,</span><span class="n">A00Real</span><span class="p">);</span></div><div class='line' id='LC424'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phip0</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">Ap0Imag</span><span class="p">,</span><span class="n">Ap0Real</span><span class="p">)</span><span class="o">-</span><span class="n">phi00</span><span class="p">;</span></div><div class='line' id='LC425'>&nbsp;&nbsp;&nbsp;<span class="n">Double_t</span> <span class="n">phim0</span><span class="o">=</span><span class="n">atan2</span><span class="p">(</span><span class="n">Am0Imag</span><span class="p">,</span><span class="n">Am0Real</span><span class="p">)</span><span class="o">-</span><span class="n">phi00</span><span class="p">;</span></div><div class='line' id='LC426'><br/></div><div class='line' id='LC427'>&nbsp;&nbsp;&nbsp;<span class="k">switch</span><span class="p">(</span><span class="n">code</span><span class="p">)</span></div><div class='line' id='LC428'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC429'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projections to hs</span></div><div class='line' id='LC430'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">1</span>:</div><div class='line' id='LC431'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC432'><br/></div><div class='line' id='LC433'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC434'><br/></div><div class='line' id='LC435'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC436'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC437'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC438'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC439'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC440'><br/></div><div class='line' id='LC441'><br/></div><div class='line' id='LC442'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC443'><br/></div><div class='line' id='LC444'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC445'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projections to Phi1</span></div><div class='line' id='LC446'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">2</span>:</div><div class='line' id='LC447'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC448'><br/></div><div class='line' id='LC449'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC450'><br/></div><div class='line' id='LC451'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC452'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC453'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC454'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC455'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC456'><br/></div><div class='line' id='LC457'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC458'><br/></div><div class='line' id='LC459'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC460'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projections to h2</span></div><div class='line' id='LC461'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">3</span>:</div><div class='line' id='LC462'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC463'><br/></div><div class='line' id='LC464'><br/></div><div class='line' id='LC465'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC466'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC467'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="o">-</span><span class="mi">8</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC468'><br/></div><div class='line' id='LC469'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">4</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">2</span><span class="o">*</span><span class="n">h2</span><span class="o">*</span><span class="n">R2Val</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC470'><br/></div><div class='line' id='LC471'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">4</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h2</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">-</span> <span class="mi">2</span><span class="o">*</span><span class="n">h2</span><span class="o">*</span><span class="n">R2Val</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC472'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC473'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC474'><br/></div><div class='line' id='LC475'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC476'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projections to h1</span></div><div class='line' id='LC477'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">4</span>:</div><div class='line' id='LC478'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC479'><br/></div><div class='line' id='LC480'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC481'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC482'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="o">-</span><span class="mi">8</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="p">(</span><span class="o">-</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">))</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC483'><br/></div><div class='line' id='LC484'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">4</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">-</span> <span class="mi">2</span><span class="o">*</span><span class="n">h1</span><span class="o">*</span><span class="n">R1Val</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC485'><br/></div><div class='line' id='LC486'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">4</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mi">1</span> <span class="o">+</span> <span class="n">Power</span><span class="p">(</span><span class="n">h1</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">2</span><span class="o">*</span><span class="n">h1</span><span class="o">*</span><span class="n">R1Val</span><span class="p">))</span><span class="o">/</span><span class="mf">3.</span><span class="p">;</span></div><div class='line' id='LC487'><br/></div><div class='line' id='LC488'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC489'><br/></div><div class='line' id='LC490'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC491'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projections to Phi</span></div><div class='line' id='LC492'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">5</span>:</div><div class='line' id='LC493'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC494'><br/></div><div class='line' id='LC495'><br/></div><div class='line' id='LC496'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC497'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC498'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC499'><br/></div><div class='line' id='LC500'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC501'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC502'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">16</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Pi</span><span class="p">())</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC503'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC504'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">f00</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fp0</span><span class="p">)</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">3</span><span class="p">)</span><span class="o">*</span><span class="n">R1Val</span><span class="o">*</span><span class="n">R2Val</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="n">Phi</span> <span class="o">+</span> <span class="n">phip0</span><span class="p">))</span><span class="o">/</span><span class="mf">4.</span><span class="p">;</span></div><div class='line' id='LC505'><br/></div><div class='line' id='LC506'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">f00</span><span class="p">)</span><span class="o">*</span><span class="n">Sqrt</span><span class="p">(</span><span class="n">fm0</span><span class="p">)</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">3</span><span class="p">)</span><span class="o">*</span><span class="n">R1Val</span><span class="o">*</span><span class="n">R2Val</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="n">Phi</span> <span class="o">-</span> <span class="n">phim0</span><span class="p">))</span><span class="o">/</span><span class="mf">4.</span><span class="p">;</span></div><div class='line' id='LC507'><br/></div><div class='line' id='LC508'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">8</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">fm0</span><span class="p">)</span><span class="o">*</span><span class="n">sqrt</span><span class="p">(</span><span class="n">fp0</span><span class="p">)</span><span class="o">*</span><span class="n">Pi</span><span class="p">()</span><span class="o">*</span><span class="n">Cos</span><span class="p">(</span><span class="mi">2</span><span class="o">*</span><span class="n">Phi</span> <span class="o">-</span> <span class="n">phim0</span> <span class="o">+</span> <span class="n">phip0</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC509'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC510'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC511'><br/></div><div class='line' id='LC512'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC513'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// projected everything</span></div><div class='line' id='LC514'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">case</span> <span class="mi">6</span>:</div><div class='line' id='LC515'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC516'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="kt">double</span> <span class="n">value</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC517'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC518'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">32</span><span class="o">*</span><span class="n">f00</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC519'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC520'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">32</span><span class="o">*</span><span class="n">fp0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC521'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC522'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="n">value</span> <span class="o">+=</span> <span class="p">(</span><span class="mi">32</span><span class="o">*</span><span class="n">fm0</span><span class="o">*</span><span class="n">Power</span><span class="p">(</span><span class="n">Pi</span><span class="p">(),</span><span class="mi">2</span><span class="p">))</span><span class="o">/</span><span class="mf">9.</span><span class="p">;</span></div><div class='line' id='LC523'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </div><div class='line' id='LC524'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span class="k">return</span> <span class="n">value</span><span class="p">;</span></div><div class='line' id='LC525'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC526'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC527'>&nbsp;&nbsp;&nbsp;<span class="n">assert</span><span class="p">(</span><span class="mi">0</span><span class="p">)</span> <span class="p">;</span></div><div class='line' id='LC528'>&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="mi">0</span> <span class="p">;</span></div><div class='line' id='LC529'><span class="p">}</span></div><div class='line' id='LC530'><br/></div><div class='line' id='LC531'><br/></div><div class='line' id='LC532'><br/></div></pre></div>
            </td>
          </tr>
        </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.14813s from github-fe137-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/awhitbeck/usercode/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

  </body>
</html>

