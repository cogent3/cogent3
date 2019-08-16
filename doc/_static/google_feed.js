google.load("feeds", "1");

function initialize() {
  var feed = new google.feeds.Feed("http://pycogent.wordpress.com/feed/");
  feed.load(function(result) {
    if (!result.error) {
      var container = document.getElementById("feed");
      for (var i = 0; i < 3; i++) {
        var entry = result.feed.entries[i];
        var tr =document.createElement('tr');
        var td =document.createElement('td');
        var link = document.createElement('a');
        link.setAttribute('href', entry.link);
        var str = entry.publishedDate;
        var patt1 = /[0-9]{2} \w+ [0-9]{4}/i;
        var pubdate=str.match(patt1).toString();
        var splitdate=pubdate.split(" ");
        var title=document.createTextNode(entry.title);
        var subtitle = document.createElement('b');
        subtitletxt=document.createTextNode(' ('+splitdate[1]+' '+splitdate[0]+', '+splitdate[2]+')');
        subtitle.setAttribute('style','color:white;font-size:9px;') 
        subtitle.appendChild(subtitletxt);
        link.appendChild(title);
        link.appendChild(subtitle)
        td.appendChild(link);
        tr.appendChild(td);
        container.appendChild(tr);
      }
    }
  });
}
google.setOnLoadCallback(initialize);