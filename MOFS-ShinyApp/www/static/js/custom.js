jQuery(function($){'use strict';$(window).on('scroll',function(){if($(this).scrollTop()>50){$('.main-nav').addClass('menu-shrink');}else{$('.main-nav').removeClass('menu-shrink');}});jQuery('.mean-menu').meanmenu({meanScreenWidth:"991"});$('#Container').mixItUp();$('.engineer-slider').owlCarousel({items:4,loop:true,margin:15,singleItem:true,nav:false,dots:true,smartSpeed:1000,autoplay:true,autoplayTimeout:4000,autoplayHoverPause:true,responsive:{0:{items:1,},600:{items:2,},1000:{items:4,}}});$('.review-slider').owlCarousel({items:1,loop:true,margin:15,nav:true,dots:false,smartSpeed:1000,autoplay:true,autoplayTimeout:4000,autoplayHoverPause:true,navText:["<i class='bx bx-chevron-left'></i>","<i class='bx bx-chevron-right'></i>"],});$('.popup-youtube').magnificPopup({disableOn:320,type:'iframe',mainClass:'mfp-fade',removalDelay:160,preloader:false,fixedContentPos:false});$('.banner-slider').owlCarousel({items:1,loop:true,margin:15,nav:true,dots:true,smartSpeed:1000,autoplay:true,autoplayTimeout:4000,autoplayHoverPause:true,navText:["<i class='bx bx-chevron-left'></i>","<i class='bx bx-chevron-right'></i>"],});$('.service-details-slider').owlCarousel({items:1,loop:true,margin:15,nav:false,dots:true,smartSpeed:1000,autoplay:true,autoplayTimeout:4000,autoplayHoverPause:true,});new WOW().init();$(".newsletter-form").validator().on("submit",function(event){if(event.isDefaultPrevented()){formErrorSub();submitMSGSub(false,"Please enter your email correctly.");}else{event.preventDefault();}});function callbackFunction(resp){if(resp.result==="success"){formSuccessSub();}
else{formErrorSub();}}
function formSuccessSub(){$(".newsletter-form")[0].reset();submitMSGSub(true,"Thank you for subscribing!");setTimeout(function(){$("#validator-newsletter").addClass('hide');},4000)}
function formErrorSub(){$(".newsletter-form").addClass("animated shake");setTimeout(function(){$(".newsletter-form").removeClass("animated shake");},1000)}
function submitMSGSub(valid,msg){if(valid){var msgClasses="validation-success";}else{var msgClasses="validation-danger";}
$("#validator-newsletter").removeClass().addClass(msgClasses).text(msg);}
$(".newsletter-form").ajaxChimp({url:"https://envytheme.us20.list-manage.com/subscribe/post?u=60e1ffe2e8a68ce1204cd39a5&amp;id=42d6d188d9",callback:callbackFunction});$('.accordion > li:eq(0) a').addClass('active').next().slideDown();$('.accordion a').on('click',function(j){var dropDown=$(this).closest('li').find('p');$(this).closest('.accordion').find('p').not(dropDown).slideUp();if($(this).hasClass('active')){$(this).removeClass('active');}else{$(this).closest('.accordion').find('a.active').removeClass('active');$(this).addClass('active');}
dropDown.stop(false,true).slideToggle();j.preventDefault();});let getDaysId=document.getElementById('days');if(getDaysId!==null){const second=1000;const minute=second*60;const hour=minute*60;const day=hour*24;let countDown=new Date('July 30, 2021 00:00:00').getTime();setInterval(function(){let now=new Date().getTime();let distance=countDown-now;document.getElementById('days').innerText=Math.floor(distance/(day)),document.getElementById('hours').innerText=Math.floor((distance%(day))/(hour)),document.getElementById('minutes').innerText=Math.floor((distance%(hour))/(minute)),document.getElementById('seconds').innerText=Math.floor((distance%(minute))/second);},second);};jQuery(window).on('load',function(){jQuery(".loader").fadeOut(500);});$('body').append('<div id="toTop" class="back-to-top-btn"><i class="bx bx-up-arrow"></i></div>');$(window).scroll(function(){if($(this).scrollTop()!=0){$('#toTop').fadeIn();}else{$('#toTop').fadeOut();}});$('#toTop').on('click',function(){$("html, body").animate({scrollTop:0},900);return false;});}(jQuery));